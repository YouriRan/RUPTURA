#!/usr/bin/env python3
"""Apply the complete-paper reference-data update to create.ipynb."""

import json
from pathlib import Path


path = Path(__file__).with_name("create.ipynb")
notebook = json.loads(path.read_text())


def set_source(index, source):
    notebook["cells"][index]["source"] = source.splitlines(keepends=True)


set_source(
    0,
    """# Marx binary and ternary breakthrough/temperature comparisons

This single plotting notebook shows Figure 16 and Figure 17 first, followed by every binary and ternary simulation in consecutive six-run batches. Ruptura results use solid lines and every native output point over the complete simulated duration; all 60 PDF-digitized references use their original 50 open-marker points. No curve is interpolated. All dimensions use SI units, and every figure is saved as PDF and PNG at 12 in width.
""",
)


imports = "".join(notebook["cells"][2]["source"])
imports = imports.replace(
    'REFERENCE_ROOT = MARX2_DIR.parents[2] / "marx" / "hassan_data"',
    'REFERENCE_ROOT = MARX2_DIR / "reference_digitized"',
)
imports = imports.replace(
    'FIGURE_DIR = MARX2_DIR / "figures"\nFIGURE_DIR.mkdir(exist_ok=True)',
    'PDF_DIR = MARX2_DIR / "output" / "pdf"\n'
    'PNG_DIR = MARX2_DIR / "figures"\n'
    'PDF_DIR.mkdir(parents=True, exist_ok=True)\n'
    'PNG_DIR.mkdir(exist_ok=True)',
)
imports = imports.replace('print(f"Figure output:   {FIGURE_DIR}")', 'print(f"PDF output:      {PDF_DIR}")\nprint(f"PNG output:      {PNG_DIR}")')
imports = imports.replace("N_PLOT_POINTS = 50", "N_REFERENCE_POINTS = 50")
set_source(2, imports)

registry = r'''REFERENCE_RUNS = {
    **{f"B{index}": {"reference": f"B{index}"} for index in range(1, 37)},
    **{f"T{index}": {"reference": f"T{index}"} for index in range(1, 25)},
}

FIGURE_16 = ["T15", "T3", "T13", "T21"]
FIGURE_17 = ["T17", "T5", "T23", "T11"]
BINARY_BATCHES = [
    [f"B{index}" for index in range(start, start + 6)]
    for start in range(1, 37, 6)
]
TERNARY_BATCHES = [
    [f"T{index}" for index in range(start, start + 6)]
    for start in range(1, 25, 6)
]


def inferred_flow_rate_m3_s(velocity_m_s):
    candidates = {1.0e-5: 0.027451748447320297, 2.0e-5: 0.054903496894640594, 3.0e-5: 0.08235524534196088}
    return min(candidates, key=lambda flow: abs(candidates[flow] - velocity_m_s))


def run_metadata(run_id):
    kind = "binary" if run_id.startswith("B") else "ternary"
    sim_dir = MARX2_DIR / kind / run_id
    with (sim_dir / "simulation.json").open() as handle:
        settings = json.load(handle)

    components = []
    for component_index, component in enumerate(settings["Components"]):
        if component.get("CarrierGas", False):
            continue
        components.append({
            "index": component_index,
            "name": component["Name"],
            "fraction_percent": float(component["GasPhaseMolFraction"]),
            "filename": f"component_{component_index}_{component['Name']}.data",
        })

    flow_m3_s = inferred_flow_rate_m3_s(float(settings["ColumnEntranceVelocity"]))
    metadata = {
        "kind": kind,
        "sim_dir": sim_dir,
        "components": components,
        "temperature_K": float(settings["InfluxTemperature"]),
        "flow_m3_s": flow_m3_s,
        "pressure_Pa": float(settings["InletPressure"]),
        "simulation_end_s": float(settings["SimulationTime"]),
        "expected_output_points": int(float(settings["SimulationTime"]) // float(settings["WriteInterval"])),
    }
    metadata.update(REFERENCE_RUNS[run_id])
    return metadata


RUNS = {
    run_id: run_metadata(run_id)
    for run_id in [
        *[f"B{index}" for index in range(1, 37)],
        *[f"T{index}" for index in range(1, 25)],
    ]
}

assert len(RUNS) == 60
assert len(BINARY_BATCHES) == 6
assert len(TERNARY_BATCHES) == 4
print(f"Registered {len(RUNS)} simulations and {len(REFERENCE_RUNS)} reference datasets.")
'''
set_source(4, registry)

readers = r'''@lru_cache(maxsize=None)
def data_column_index(path_string, heading):
    """Resolve a one-based Ruptura '# column N:' header to a NumPy index."""
    prefix = heading.casefold()
    with Path(path_string).open("rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.startswith("# column "):
                number_text, description = line[9:].split(":", 1)
                if description.strip().casefold().startswith(prefix):
                    return int(number_text) - 1
    raise ValueError(f"Column heading {heading!r} not found in {path_string}")


@lru_cache(maxsize=None)
def outlet_rows(path_string, min_columns):
    # Parse only the final axial row before each blank-line block separator.
    # The memory map keeps even multi-gigabyte histories memory-efficient.
    rows = []
    with Path(path_string).open("rb") as handle:
        with mmap.mmap(handle.fileno(), length=0, access=mmap.ACCESS_READ) as mapped:
            block_start = 0
            while True:
                block_end = mapped.find(b"\n\n", block_start)
                if block_end < 0:
                    break
                line_start = mapped.rfind(b"\n", block_start, block_end) + 1
                line = mapped[line_start:block_end].strip()
                if line and not line.startswith(b"#"):
                    values = np.fromstring(line.decode("ascii"), sep=" ")
                    if len(values) >= min_columns:
                        rows.append(values[:min_columns])
                block_start = block_end + 2
    if not rows:
        raise ValueError(f"No complete {min_columns}-column blocks found in {path_string}")
    return np.asarray(rows, dtype=float)


def direct_xy(x, y, source):
    """Return source samples unchanged after basic integrity checks."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError(f"{source}: x and y must be matching one-dimensional arrays")
    if len(x) < 2 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError(f"{source}: curve is incomplete or contains non-finite values")
    if np.any(np.diff(x) < 0.0):
        raise ValueError(f"{source}: time values are not ordered")
    return x, y


def load_reference(run_id):
    settings = RUNS[run_id]
    ref_dir = REFERENCE_ROOT / settings["reference"]
    reference = {"temperature": {}}
    for component in settings["components"]:
        species = component["name"]
        raw = np.loadtxt(ref_dir / f"exp-{species}.csv", delimiter=",")
        reference[species] = direct_xy(raw[:, 0], raw[:, 1], f"{run_id} reference {species}")
    raw_temperature = np.loadtxt(ref_dir / "exp-Temp.csv", delimiter=",")
    reference["temperature"]["gas"] = direct_xy(
        raw_temperature[:, 0], raw_temperature[:, 1], f"{run_id} reference temperature"
    )
    return reference


@lru_cache(maxsize=None)
def load_run(run_id):
    settings = RUNS[run_id]
    reference = load_reference(run_id)

    sim_dir = settings["sim_dir"]
    column_path = (sim_dir / "column.data").resolve()
    column_time_index = data_column_index(str(column_path), "Time, t")
    temperature_indices = {
        "gas": data_column_index(str(column_path), "Gas temperature,"),
        "solid": data_column_index(str(column_path), "Solid temperature,"),
        "wall": data_column_index(str(column_path), "Wall temperature,"),
    }
    column = outlet_rows(str(column_path), min_columns=max(temperature_indices.values()) + 1)
    time_s = 60.0 * column[:, column_time_index]
    if len(time_s) != settings["expected_output_points"]:
        raise RuntimeError(
            f"{run_id}: expected {settings['expected_output_points']} column outputs, "
            f"found {len(time_s)}. Rerun this case."
        )
    time_s, _ = direct_xy(time_s, time_s, f"{run_id} column time")

    simulation = {"temperature": {}}
    for component in settings["components"]:
        component_path = (sim_dir / component["filename"]).resolve()
        component_time_index = data_column_index(str(component_path), "Time, t")
        mole_fraction_index = data_column_index(str(component_path), "Mole fraction,")
        rows = outlet_rows(str(component_path), min_columns=max(component_time_index, mole_fraction_index) + 1)
        component_time_s = 60.0 * rows[:, component_time_index]
        if len(component_time_s) != settings["expected_output_points"]:
            raise RuntimeError(
                f"{run_id}/{component['name']}: expected {settings['expected_output_points']} outputs, "
                f"found {len(component_time_s)}."
            )
        simulation[component["name"]] = direct_xy(
            component_time_s,
            rows[:, mole_fraction_index],
            f"{run_id} simulation {component['name']}",
        )
    for label, column_index in temperature_indices.items():
        simulation["temperature"][label] = direct_xy(
            time_s, column[:, column_index], f"{run_id} simulation {label} temperature"
        )
    return simulation, reference


for run_id in RUNS:
    simulation, reference = load_run(run_id)
    for component in RUNS[run_id]["components"]:
        species = component["name"]
        assert len(simulation[species][0]) == RUNS[run_id]["expected_output_points"]
        assert len(reference[species][0]) == N_REFERENCE_POINTS
    for phase in ("gas", "solid", "wall"):
        assert len(simulation["temperature"][phase][0]) == RUNS[run_id]["expected_output_points"]
    assert len(reference["temperature"]["gas"][0]) == N_REFERENCE_POINTS

print(
    f"Validated {len(RUNS)} simulations using every native Ruptura output point; "
    f"each reference curve retains its {N_REFERENCE_POINTS} digitized points."
)
'''
set_source(6, readers)

plotting = "".join(notebook["cells"][8]["source"])
plotting = plotting.replace(
    '"N2":  {"color": "C1", "marker": "o", "label": r"N$_2$"},\n    "CO2": {"color": "C2", "marker": "D", "label": r"CO$_2$"},',
    '"N2":  {"color": "C1", "marker": "D", "label": r"N$_2$"},\n    "CO2": {"color": "C2", "marker": "o", "label": r"CO$_2$"},',
)
plotting = plotting.replace('include_reference = any("reference" in RUNS[run_id] for run_id in run_ids)', 'include_reference = True')
plotting = plotting.replace('pdf_path = FIGURE_DIR / f"{stem}.pdf"', 'pdf_path = PDF_DIR / f"{stem}.pdf"')
plotting = plotting.replace('png_path = FIGURE_DIR / f"{stem}.png"', 'png_path = PNG_DIR / f"{stem}.png"')
plotting = plotting.replace(
    'fig = plt.figure(figsize=(FIGURE_WIDTH, 8.7), constrained_layout=True)',
    'figure_height = 8.7 if ncols == 2 else 10.0\n    fig = plt.figure(figsize=(FIGURE_WIDTH, figure_height), constrained_layout=True)',
)
plotting = plotting.replace(
    'width_ratios=[1.0] * ncols + [0.42], wspace=0.12, hspace=0.12',
    'width_ratios=[1.0] * ncols + [0.62], wspace=0.10, hspace=0.12',
)
set_source(8, plotting)

validation = '''expected_stems = ['figure_16', 'figure_17', 'B1_B6', 'B7_B12', 'B13_B18', 'B19_B24', 'B25_B30', 'B31_B36', 'T1_T6', 'T7_T12', 'T13_T18', 'T19_T24']

for stem in expected_stems:
    for directory, suffix in ((PDF_DIR, '.pdf'), (PNG_DIR, '.png')):
        path = directory / f'{stem}{suffix}'
        assert path.is_file() and path.stat().st_size > 0, path
        print(f'{path.name:28s} {path.stat().st_size / 1024:9.1f} KiB')
'''
set_source(34, validation)

notebook["cells"][3]["source"] = [
    "## 2. Run registry and figure groups\n",
    "\n",
    "All 60 runs use experimental reference traces digitized from the paper PDFs. Figure 16 and Figure 17 remain first; subsequent sections group B and T runs in batches of six.\n",
]
notebook["cells"][5]["source"] = [
    "## 3. Native output readers (no interpolation)\n",
    "\n",
    "Every complete Ruptura outlet history is read and plotted at its native write frequency and full simulated duration. The PDF-digitized references retain their original 50 points.\n",
]

for cell in notebook["cells"]:
    if cell["cell_type"] == "code":
        cell["execution_count"] = None
        cell["outputs"] = []

path.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Updated {path}")
