"""Gas/liquid FDM/FVM breakthrough and cycles using the shared native column model."""
from __future__ import annotations
import copy
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

def _unique_pairs(pairs):
    result={}
    for key,value in pairs:
        if key in result:raise ValueError(f"Duplicate JSON key: {key}")
        result[key]=value
    return result

def canonical_chemisorption(document):
    """Import the retired chemistry block using the native canonical converter."""
    def legacy(value):
        if isinstance(value, dict):
            return any(k.lower() == "chemicalkinetics" or legacy(v) for k, v in value.items())
        return isinstance(value, list) and any(legacy(v) for v in value)
    if not legacy(document):
        return copy.deepcopy(document)
    try:
        from ._ruptura import canonicalize_chemisorption_json
    except ImportError:
        candidate = os.environ.get("RUPTURA_CLI_EXECUTABLE") or shutil.which("ruptura") or str(Path(__file__).resolve().parents[1]/"build/ruptura")
        if not Path(candidate).is_file():
            raise ValueError("Build the current RUPTURA engine to convert legacy ChemicalKinetics into ChemisorptionSites.")
        with tempfile.TemporaryDirectory(prefix="ruptura-chemistry-import-") as directory:
            path=Path(directory)/"simulation.json";path.write_text(json.dumps(document,allow_nan=False))
            run=subprocess.run([candidate,str(path),"--canonicalize"],capture_output=True,text=True,timeout=30)
            if run.returncode: raise ValueError(run.stderr or run.stdout)
            return json.loads(run.stdout)
    return json.loads(canonicalize_chemisorption_json(json.dumps(document,allow_nan=False)))


class ProcessSimulation:
    """Configure and run a main-style input without mutating the source file.

    The same Components, Geometry and material definitions apply to breakthrough
    and cycles. Direct calls execute the C++ solver; Python does not integrate.
    """
    def __init__(self, input, directory=None):
        if isinstance(input, (str, Path)):
            path = Path(input).resolve()
            self.input = json.loads(path.read_text(),object_pairs_hook=_unique_pairs)
            self.directory = Path(directory).resolve() if directory else path.parent
        elif isinstance(input, dict):
            self.input = copy.deepcopy(input)
            self.directory = Path(directory or ".").resolve()
        else:
            raise TypeError("input must be a simulation.json path or a dictionary")
        self.input = canonical_chemisorption(self.input)
        self.input.setdefault("SpatialDiscretization", "FVM")

    def validate(self):
        from ._ruptura import validate_process_json
        return json.loads(validate_process_json(json.dumps(self.input), str(self.directory)))

    def run(self, *, restart=False, cycles=None, plots=None):
        from ._ruptura import run_process_json
        config = copy.deepcopy(self.input)
        if restart:
            config["Restart"] = True
        if cycles is not None:
            if not isinstance(cycles, int) or cycles < 1:
                raise ValueError("cycles must be a positive integer")
            config["NumberOfCycles"] = cycles
            config.setdefault("CyclicSteadyState", {})["Enabled"] = False
        if plots is not None:
            config.setdefault("Output", {})["Plots"] = bool(plots)
        return json.loads(run_process_json(json.dumps(config), str(self.directory)))

    @staticmethod
    def capabilities():
        from ._ruptura import process_capabilities_json
        return json.loads(process_capabilities_json())

    def results(self):
        """Return saved NumPy fields, the summary, and their original metadata."""
        import numpy as np
        output = Path(self.input.get("Output", {}).get("Directory", "results"))
        if not output.is_absolute():
            output = self.directory / output
        with np.load(output / "states.npz", allow_pickle=False) as arrays:
            result = {key: arrays[key].copy() for key in arrays.files}
        result["summary"] = json.loads((output / "summary.json").read_text())
        result["metadata"] = json.loads((output / "metadata.json").read_text())
        return result

    def compute(self):
        self.run(plots=False)
        return self.results()["state"]
