"""Persistent local jobs. Only owned CLI processes are ever canceled."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import copy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import threading
import uuid

from .config import PAPER_PRESETS, import_document, migrate_example_catalog, normalize_input, presets


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)


class StateStore:
    def __init__(self, workspace):
        self.path = Path(workspace) / ".ruptura-lab" / "state.json"
        self.lock = threading.RLock()
        migrated = False
        if self.path.exists():
            loaded = json.loads(self.path.read_text())
            self.state = import_document(loaded)
            migrated = loaded != self.state
        else:
            identifier = uuid.uuid4().hex
            self.state = {"Format": "RupturaWebLabState", "Version": 1, "selectedId": identifier,
                          "simulations": [{"id": identifier, "input": presets()["breakthrough"]}]}
        # Seed new paper examples once, preserving edits, selection, and deliberate deletions.
        included = set(self.state.get("includedPresets", []))
        catalog = presets()
        changed = migrated or not self.path.exists()
        for key in PAPER_PRESETS:
            if key in included:
                continue
            if not any(item.get("presetId") == key or item["id"] == f"default-{key}" for item in self.state["simulations"]):
                self.state["simulations"].append({"id": f"default-{key}", "presetId": key, "input": catalog[key]})
            included.add(key)
            changed = True
        changed = migrate_example_catalog(self.state) or changed
        before = self.state
        self.state = import_document(self.state)
        changed = changed or self.state != before
        if changed:
            self.state["includedPresets"] = sorted(included)
            atomic_json(self.path, self.state)

    def get(self):
        with self.lock:
            return copy.deepcopy(self.state)

    def save(self, state):
        value = import_document(state)
        with self.lock:
            atomic_json(self.path, value)
            self.state = value
        return self.get()


def find_cli(explicit=None):
    configured = explicit or os.environ.get("RUPTURA_CLI_EXECUTABLE")
    if configured:
        executable = shutil.which(str(configured))
        if not executable:
            raise ValueError(f"Ruptura executable not found or not executable: {configured}")
        return str(Path(executable).resolve())
    for candidate in (shutil.which("ruptura"), Path(__file__).resolve().parents[2] / "build" / "ruptura",
                      Path(__file__).resolve().parents[2] / "build" / "conda" / "ruptura"):
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return str(Path(candidate).resolve())
    return None


class RunManager:
    def __init__(self, workspace, cli=None):
        self.workspace = Path(workspace).resolve()
        self.root = self.workspace / "simulations"
        self.root.mkdir(parents=True, exist_ok=True)
        self.cli = find_cli(cli)
        self.lock = threading.RLock()
        self.records = {}
        self.processes = {}
        self.futures = {}
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="ruptura-run")
        self.closed = False
        for path in self.root.glob("*/manifest.json"):
            try:
                record = json.loads(path.read_text())
                identifier = path.parent.name
                record.update(id=identifier, directory=str(path.parent.resolve()), simulationPath=str(path.parent / "simulation.json"),
                              notebookPath=str(path.parent / "analysis.ipynb"), logPath=str(path.parent / "run.log"))
                if record.get("status") in ("running", "queued", "canceling"):
                    record.update(status="interrupted", errorMessage="The Lab server stopped before this run completed.")
                    atomic_json(path, record)
                self.records[identifier] = record
            except (OSError, ValueError, TypeError):
                continue

    def capabilities(self):
        result = {"cli": self.cli, "sundials": False, "input_validation": False}
        if not self.cli:
            result["error"] = "Build ruptura or choose its executable in Settings."
            return result
        try:
            process = subprocess.run([self.cli, "--capabilities"], capture_output=True, text=True, timeout=10)
            if process.returncode:
                raise ValueError("The selected executable needs rebuilding with --capabilities and --validate support.")
            result.update(json.loads(process.stdout))
        except (OSError, ValueError, subprocess.TimeoutExpired) as error:
            result["error"] = str(error)
        return result

    def validate(self, document, directory=None):
        config = normalize_input(document, directory or self.workspace)
        if not self.cli:
            raise ValueError("Ruptura was not found. Choose the executable in Settings, or build the CLI first.")
        capabilities = self.capabilities()
        if not capabilities.get("input_validation"):
            raise ValueError(capabilities.get("error", "The CLI does not support validation."))
        if config.get("BreakthroughIntegrator") == "CVODE" and not capabilities.get("sundials"):
            raise ValueError("This CLI was built without CVODE. Select RungeKutta3 or SIRK3.")
        with tempfile.TemporaryDirectory(prefix="ruptura-validate-") as temporary:
            path = Path(temporary) / "simulation.json"
            path.write_text(json.dumps(config, allow_nan=False))
            process = subprocess.run([self.cli, "--validate", str(path)], cwd=self.workspace, capture_output=True, text=True, timeout=30)
            if process.returncode:
                raise ValueError((process.stderr or process.stdout or "Invalid simulation input.").strip())
        return config

    def list(self):
        with self.lock:
            return sorted((copy.deepcopy(item) for item in self.records.values()), key=lambda item: item["id"], reverse=True)

    def get(self, identifier):
        with self.lock:
            if identifier not in self.records:
                raise FileNotFoundError("Run not found.")
            return copy.deepcopy(self.records[identifier])

    def persist(self, record):
        atomic_json(Path(record["directory"]) / "manifest.json", record)

    def start(self, document, simulation_id=None):
        config = self.validate(document)
        name = config.get("DisplayName", "simulation")
        slug = re.sub(r"[^A-Za-z0-9._-]+", "-", str(name)).strip("-.")[:32] or "simulation"
        identifier = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8] + "-" + slug
        directory = self.root / identifier
        with self.lock:
            if self.closed:
                raise ValueError("The Lab server is stopping.")
            directory.mkdir()
            atomic_json(directory / "simulation.json", config)
            (directory / "run.log").write_text(f"Ruptura: {self.cli}\nWorking directory: {directory}\n\n")
            record = {"id": identifier, "name": name, "simulationId": simulation_id, "simulationType": config.get("SimulationType", "Breakthrough"),
                      "directory": str(directory), "simulationPath": str(directory / "simulation.json"), "notebookPath": str(directory / "analysis.ipynb"),
                      "logPath": str(directory / "run.log"), "status": "queued", "createdAt": timestamp(), "startedAt": None,
                      "finishedAt": None, "exitCode": None, "errorMessage": ""}
            self.records[identifier] = record
            self.persist(record)
            self.futures[identifier] = self.executor.submit(self._execute, identifier, self.cli)
        return self.get(identifier)

    def _execute(self, identifier, executable):
        process = None
        try:
            with self.lock:
                record = self.records[identifier]
                if record["status"] == "canceled":
                    return
                log = open(record["logPath"], "ab", buffering=0)
                try:
                    process = subprocess.Popen([executable], cwd=record["directory"], stdin=subprocess.DEVNULL,
                                               stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                finally:
                    log.close()
                self.processes[identifier] = process
                record.update(status="running", startedAt=timestamp())
                self.persist(record)
            code = process.wait()
            with self.lock:
                record = self.records[identifier]
                canceled = record["status"] in ("canceling", "canceled")
                record.update(status="canceled" if canceled else "completed" if code == 0 else "failed", finishedAt=timestamp(), exitCode=code)
                if not canceled and code != 0:
                    record["errorMessage"] = self.log(identifier)["text"][-2000:] or f"ruptura exited with code {code}."
                self.persist(record)
        except Exception as error:
            with self.lock:
                record = self.records[identifier]
                record.update(status="failed", finishedAt=timestamp(), errorMessage=str(error))
                self.persist(record)
        finally:
            with self.lock:
                self.processes.pop(identifier, None)

    def cancel(self, identifier):
        with self.lock:
            record = self.records.get(identifier)
            if not record:
                raise FileNotFoundError("Run not found.")
            if record["status"] not in ("queued", "running", "canceling"):
                raise ValueError("This run has already finished.")
            process = self.processes.get(identifier)
            if process is None:
                future = self.futures.get(identifier)
                if future:
                    future.cancel()
                record.update(status="canceled", finishedAt=timestamp(), errorMessage="Canceled by user.")
            else:
                record["status"] = "canceling"
                if process.poll() is None:
                    process.terminate()
                timer = threading.Timer(1.5, lambda: process.kill() if process.poll() is None else None)
                timer.daemon = True
                timer.start()
            self.persist(record)
        return self.get(identifier)

    def log(self, identifier, offset=None):
        record = self.get(identifier)
        path = Path(record["logPath"])
        if not path.exists():
            return {"text": "", "offset": 0}
        with path.open("rb") as handle:
            size = path.stat().st_size
            position = max(0, size - 64000) if offset is None else max(0, min(int(offset), size))
            handle.seek(position)
            data = handle.read(64000)
            return {"text": data.decode("utf-8", errors="replace"), "offset": handle.tell()}

    def run_file(self, identifier, relative):
        root = Path(self.get(identifier)["directory"]).resolve()
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise FileNotFoundError("Output file not found.")
        return path

    def close(self):
        with self.lock:
            self.closed = True
            active = [item["id"] for item in self.records.values() if item.get("status") in ("queued", "running", "canceling")]
        for identifier in active:
            self.cancel(identifier)
        self.executor.shutdown(wait=True, cancel_futures=True)
