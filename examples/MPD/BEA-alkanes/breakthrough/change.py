import json
import glob

for dir in glob.glob("*_r1"):
    with open(f"{dir}/simulation.json") as f:
        settings = json.load(f)
    settings["SimulationTime"] = "auto"
    settings["WriteInterval"] = 5.0
    # CVODE chooses its own internal steps and takes no TimeStep.
    if settings.get("BreakthroughIntegrator") == "CVODE":
        settings.pop("TimeStep", None)
    else:
        settings["TimeStep"] = 5e-4
    settings["NumberOfGridPoints"] = 50
    with open(f"{dir}/simulation.json", 'w') as f:
        json.dump(settings, f, indent=4)