"""Loopback-only HTTP application; no frontend toolchain or notebook server required."""
from __future__ import annotations

import argparse
import importlib.util
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import mimetypes
from pathlib import Path
import secrets
import signal
import subprocess
import sys
import threading
from urllib.parse import parse_qs, unquote, urlparse
import webbrowser

from .about import ABOUT, example_notes
from .config import PRESET_INFO, SCHEMA, import_document, presets
from .runs import RunManager, StateStore, atomic_json, find_cli

ASSETS = Path(__file__).with_name("static")


class LabServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address, workspace, cli=None):
        self.workspace = Path(workspace).resolve()
        self.workspace.mkdir(parents=True, exist_ok=True)
        self.store = StateStore(self.workspace)
        settings_path = self.workspace / ".ruptura-lab" / "settings.json"
        saved = json.loads(settings_path.read_text()) if settings_path.exists() else {}
        self.manager = RunManager(self.workspace, cli or saved.get("cli"))
        self.session = secrets.token_urlsafe(32)
        try:
            super().__init__(address, LabHandler)
        except BaseException:
            self.manager.close()
            raise

    def server_close(self):
        self.manager.close()
        super().server_close()


class LabHandler(BaseHTTPRequestHandler):
    server_version = "RupturaLab/2.0"

    def log_message(self, *_args):
        pass

    def send(self, data, status=200, content_type="application/json; charset=utf-8", filename=None, cache=False):
        if isinstance(data, (dict, list)):
            data = json.dumps(data, allow_nan=False).encode()
        elif isinstance(data, str):
            data = data.encode()
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Cache-Control", "private, max-age=3600" if cache else "no-store")
        if filename:
            from urllib.parse import quote
            self.send_header("Content-Disposition", f"attachment; filename*=UTF-8''{quote(filename)}")
        self.end_headers()
        self.wfile.write(data)

    def guard(self, api=False):
        hosts = {f"127.0.0.1:{self.server.server_port}", f"localhost:{self.server.server_port}"}
        if self.headers.get("Host") not in hosts:
            self.send({"error": "Use the local Lab address printed by the server."}, 403)
            return False
        origin = self.headers.get("Origin")
        if origin and origin not in {"http://" + host for host in hosts}:
            self.send({"error": "Cross-origin requests are not accepted."}, 403)
            return False
        if api and not secrets.compare_digest(self.headers.get("X-Ruptura-Session", ""), self.server.session):
            self.send({"error": "Refresh the Lab page to reconnect."}, 403)
            return False
        return True

    def body(self):
        size = int(self.headers.get("Content-Length", "0"))
        if not 0 < size <= 10 * 1024 * 1024:
            raise ValueError("Request must contain JSON smaller than 10 MB.")
        result = json.loads(self.rfile.read(size), parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"Invalid number: {value}")))
        if not isinstance(result, dict):
            raise ValueError("Expected a JSON object.")
        return result

    def handle_error(self, error):
        code = 404 if isinstance(error, FileNotFoundError) else 503 if isinstance(error, subprocess.TimeoutExpired) else 400
        self.send({"error": str(error) or type(error).__name__}, code)

    def do_GET(self):
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        if not self.guard(api=path.startswith("/api/") and path != "/api/bootstrap"):
            return
        try:
            query = parse_qs(parsed.query)
            if path == "/api/bootstrap":
                missing = [name for name in ("numpy", "pandas", "plotly") if importlib.util.find_spec(name) is None]
                self.send({"session": self.server.session, "state": self.server.store.get(), "schema": SCHEMA, "presets": presets(),
                           "presetInfo": PRESET_INFO, "presetNotes": {key: example_notes({"input": value, "presetId": key}) for key, value in presets().items()},
                           "about": ABOUT, "runs": self.server.manager.list(), "capabilities": self.server.manager.capabilities(),
                           "workspace": str(self.server.workspace), "python": sys.executable, "missingDependencies": missing})
            elif path == "/api/runs":
                self.send(self.server.manager.list())
            elif path == "/api/state":
                self.send(self.server.store.get())
            elif path == "/api/files":
                directory = Path(query.get("path", [str(self.server.workspace)])[0]).expanduser().resolve()
                if not directory.is_dir():
                    raise ValueError("Choose an existing directory.")
                entries = []
                for item in sorted(directory.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower())):
                    if not item.name.startswith(".") and (item.is_dir() or item.suffix.lower() == ".json"):
                        entries.append({"name": item.name, "path": str(item), "directory": item.is_dir()})
                    if len(entries) >= 1000:
                        break
                self.send({"path": str(directory), "parent": str(directory.parent), "entries": entries})
            elif path.startswith("/api/runs/"):
                parts = path.split("/")
                identifier = parts[3]
                record = self.server.manager.get(identifier)
                if len(parts) == 4:
                    root = Path(record["directory"])
                    files = [{"name": item.name, "bytes": item.stat().st_size} for item in root.iterdir() if item.is_file() and item.name != "manifest.json"]
                    config = json.loads((root / "simulation.json").read_text())
                    self.send({"record": record, "input": config, "files": files})
                elif parts[4] == "log":
                    self.send(self.server.manager.log(identifier, query.get("offset", [None])[0]))
                elif parts[4] == "plots":
                    from .plots import metadata
                    self.send(metadata(record["directory"]))
                elif parts[4] == "file" and len(parts) == 6:
                    file = self.server.manager.run_file(identifier, parts[5])
                    self.send_file(file)
                else:
                    raise FileNotFoundError("Endpoint not found.")
            elif path == "/vendor/plotly.js":
                try:
                    from plotly.offline import get_plotlyjs
                    self.send(get_plotlyjs(), content_type="application/javascript; charset=utf-8", cache=True)
                except ImportError:
                    self.send("window.rupturaPlotlyMissing = true;", content_type="application/javascript")
            elif path in ("/", "/app.js", "/style.css", "/ruptura-logo.png"):
                file = ASSETS / ("index.html" if path == "/" else path[1:])
                self.send(file.read_bytes(), content_type=mimetypes.guess_type(file)[0] or "text/plain")
            elif path == "/favicon.ico":
                self.send(b"", 204, content_type="image/x-icon")
            else:
                raise FileNotFoundError("Page not found.")
        except (ValueError, OSError, KeyError, TypeError, ImportError, subprocess.TimeoutExpired) as error:
            self.handle_error(error)

    def send_file(self, path):
        from urllib.parse import quote
        self.send_response(200)
        self.send_header("Content-Type", mimetypes.guess_type(path)[0] or "application/octet-stream")
        self.send_header("Content-Length", str(path.stat().st_size))
        self.send_header("Content-Disposition", f"attachment; filename*=UTF-8''{quote(path.name)}")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        with path.open("rb") as handle:
            while chunk := handle.read(1 << 20):
                self.wfile.write(chunk)

    def do_POST(self):
        if not self.guard(api=True):
            return
        try:
            data = self.body()
            path = urlparse(self.path).path
            if path == "/api/state":
                self.send(self.server.store.save(data))
            elif path == "/api/import":
                if data.get("path"):
                    file = Path(data["path"]).expanduser().resolve()
                    document = json.loads(file.read_text())
                    directory = file.parent
                else:
                    document = data.get("document")
                    directory = Path(data.get("directory") or self.server.workspace).expanduser().resolve()
                if not isinstance(document, dict):
                    raise ValueError("Choose a simulation or Lab workspace JSON object.")
                imported = import_document(document, directory)
                imported["documentType"] = "workspace" if document.get("Format") or "Simulations" in document else "simulation"
                self.send(imported)
            elif path == "/api/validate":
                config = self.server.manager.validate(data.get("input"), data.get("directory"))
                self.send({"valid": True, "input": config})
            elif path == "/api/settings":
                cli = find_cli(data.get("cli"))
                previous = self.server.manager.cli
                self.server.manager.cli = cli
                capabilities = self.server.manager.capabilities()
                if capabilities.get("error"):
                    self.server.manager.cli = previous
                    raise ValueError(capabilities["error"])
                atomic_json(self.server.workspace / ".ruptura-lab" / "settings.json", {"cli": cli})
                self.send(capabilities)
            elif path == "/api/runs":
                self.send(self.server.manager.start(data.get("input"), data.get("simulationId")), 201)
            elif path.startswith("/api/runs/"):
                parts = path.split("/")
                if len(parts) != 5:
                    raise FileNotFoundError("Endpoint not found.")
                identifier, action = parts[3:]
                if action == "cancel":
                    self.send(self.server.manager.cancel(identifier))
                elif action == "plot":
                    from .plots import json_figure
                    self.send(json_figure(self.server.manager.get(identifier)["directory"], data))
                elif action == "report":
                    from .plots import report
                    self.send({"name": report(self.server.manager.get(identifier)["directory"]).name})
                elif action == "video":
                    from .video import export_mp4
                    self.send({"name": export_mp4(self.server.manager.get(identifier)["directory"], data).name})
                else:
                    raise FileNotFoundError("Endpoint not found.")
            else:
                raise FileNotFoundError("Endpoint not found.")
        except (ValueError, OSError, KeyError, TypeError, ImportError, subprocess.TimeoutExpired, IndexError) as error:
            self.handle_error(error)


def main():
    parser = argparse.ArgumentParser(description="Ruptura Lab — local simulation workspace and interactive plots")
    parser.add_argument("--workspace", type=Path, default=Path.cwd(), help="Workspace containing simulations/ and saved Lab state")
    parser.add_argument("--cli", help="Path to the ruptura executable")
    parser.add_argument("--port", type=int, default=8765, help="Local port; use 0 to choose an available port")
    parser.add_argument("--no-browser", action="store_true", help="Print the URL without opening the browser")
    args = parser.parse_args()
    try:
        server = LabServer(("127.0.0.1", args.port), args.workspace, args.cli)
    except (ValueError, OSError) as error:
        parser.exit(1, f"Could not start Ruptura Lab: {error}\n")
    url = f"http://127.0.0.1:{server.server_port}"
    def stop(_signum, _frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stop)
    print(f"Ruptura Lab: {url}\nWorkspace: {server.workspace}\nPress Ctrl+C to stop the server and cancel active runs.", flush=True)
    if not args.no_browser:
        timer = threading.Timer(0.4, lambda: webbrowser.open(url))
        timer.daemon = True
        timer.start()
    try:
        server.serve_forever(poll_interval=0.2)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
