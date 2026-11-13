#!/usr/bin/env python3
"""Headless Matplotlib figures for shared simulations and legacy breakthrough results.

The native solver writes NumPy-compatible ``states.npz``, ``phases/*.npz`` and
CSV summaries. This script only reads those files; it never integrates anything
and is not a numerical backend. It produces the same figure set as the Python
edition of the package, as PNG and PDF, plus ``plots/plot_manifest.json``.

    python3 tools/plot_cycle_results.py RESULTS_DIR [--formats png pdf] [--dpi 160]
    python3 tools/plot_cycle_results.py RESULTS_DIR --legacy-input INPUT_JSON

``RESULTS_DIR`` defaults to ``results`` in the current directory. The exit
status is 0 only when at least one figure was written, so the caller can fall
back to the built-in SVG overview when Matplotlib is unavailable.
"""
from __future__ import annotations

import argparse
from contextlib import nullcontext
import csv
import json
from pathlib import Path
import re
import sys

REQUIRED = "This script needs NumPy and Matplotlib: python3 -m pip install numpy matplotlib"


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def read_reference(output):
    """Optional plot reference values written by the native runner."""
    path = output / "metadata.json"
    if not path.is_file():
        return None
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    reference = record.get("plot_reference")
    if not isinstance(reference, dict):
        return None
    try:
        pressure = float(reference["pressure_Pa"])
        temperature = float(reference["temperature_K"])
        composition = [float(x) for x in reference["feed_composition"]]
        gas_constant = float(reference.get("R", 8.31446261815324))
    except (KeyError, TypeError, ValueError):
        return None
    if not (pressure > 0 and temperature > 0 and composition):
        return None
    doc=record.get("input_document",{})
    if doc.get("FluidPhase")=="Liquid":
        # Preserve the tuple interface while supplying the actual solute reference.
        composition=[float(x.get("LiquidPhaseConcentration",0.)) for x in doc.get("Components",[])]
        pressure=gas_constant*temperature
    return pressure, temperature, composition, gas_constant


def plot_results(output_dir, formats=("png", "pdf"), dpi=160):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SystemExit(f"{REQUIRED}\n{exc}")

    output = Path(output_dir)
    state_file = output / "states.npz"
    plots = output / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    saved = []
    report_path = output / "convergence_report.json"
    report = json.loads(report_path.read_text()) if report_path.is_file() else {}
    chosen = report.get("representative_cycle")
    plot_cycles = report.get("plot_cycles", [chosen] if chosen is not None else [])
    status = report.get("status", "CSS status unavailable")
    note = f"{status}; selected cycle {chosen}" if chosen is not None else status
    if report and not report.get("strict_css", False):
        note += " — NOT strict CSS"
    error_stats = report.get("error_statistics", {}).get("css_state_error", {})
    if error_stats.get("min") is not None:
        note += (f"\nTail CSS error min/mean/max: {error_stats['min']:.3g} / "
                 f"{error_stats['mean']:.3g} / {error_stats['max']:.3g}; pass ≤ 1")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.constrained_layout.use": True, "savefig.dpi": dpi})

    def save(fig, stem):
        fig.text(.5, .006, note, ha="center", va="bottom", fontsize=7,
                 color="black" if report.get("strict_css") or metadata.get("input_document",{}).get("SimulationType")=="Breakthrough" else "#8a341c")
        footer = min(.25, (.18 + .12 * len(note.splitlines())) / fig.get_figheight())
        fig.get_layout_engine().set(rect=(0, footer, 1, 1-footer), h_pad=.06, w_pad=.06)
        for suffix in formats:
            path = plots / f"{stem}.{suffix}"
            fig.savefig(path)
            saved.append(str(path))
        plt.close(fig)

    reference_values = read_reference(output)
    phase_records = []
    names = []
    trend_path = output / "history_trends.csv"
    trends = np.genfromtxt(trend_path, delimiter=",", names=True, ndmin=1) if trend_path.is_file() and trend_path.stat().st_size else None
    metadata_path = output / "metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    if metadata.get("input_document",{}).get("SimulationType")=="Breakthrough":note="Breakthrough transient; cyclic steady state is not applicable. Driving concentrations use the stated feed reference."
    liquid=metadata.get("input_document",{}).get("FluidPhase")=="Liquid"
    if liquid:note+="\nLiquid concentrations; any purity metric is a tracked-solute fraction (implicit solvent excluded)."
    if metadata.get("output_recording", {}).get("profile") == "compact":
        note += "\nSampled profiles; flux totals use solver-integrated ledgers."
    if metadata.get("output_recording", {}).get("profile_history_incomplete"):
        note += "\nINCOMPLETE profile history: some earlier profiles are unavailable."
    cycle_starts = {}
    phase_csv = output / "phase_metrics.csv"
    if phase_csv.is_file():
        with phase_csv.open(newline="") as stream:
            for row in csv.DictReader(stream):
                cyc = int(row["cycle"])
                cycle_starts[cyc] = min(cycle_starts.get(cyc, float("inf")), float(row["start_time_s"]))
    if state_file.is_file():
        with np.load(state_file, allow_pickle=False) as data:
            t = np.asarray(data["time" if "time" in data else "globaltime"], dtype=float)
            names = [str(x) for x in data["column_names"]]
            components = [str(x) for x in data["comp_names"]]
            cycle_ids = np.asarray(data["cycle"], dtype=int)
            selected_indices = np.flatnonzero(cycle_ids == chosen) if chosen is not None else np.array([])
            if chosen is not None and not selected_indices.size:
                raise ValueError(f"Selected cycle {chosen} is absent from the profile archive; it cannot be plotted as a substitute cycle.")
            selected_end = int(selected_indices[-1]) if selected_indices.size else -1
            ht = trends["time"] if trends is not None else t
            def trend(name, field, fallback):
                key = f"{name}_{field}"
                return trends[key] if trends is not None and key in trends.dtype.names else fallback
            # End-of-phase plots are views of the canonical archive, not duplicate files.
            for cyc in (plot_cycles or [int(cycle_ids[-1])]):
                for ph in np.unique(data["phase"][cycle_ids == cyc]):
                    idx = np.flatnonzero((cycle_ids == cyc) & (data["phase"] == ph))
                    if not idx.size: continue
                    record = {}
                    for name in names:
                        for field in ("c", "q", "P", "Tg", "Ts", "Tw", "Tj"):
                            key = f"{name}_{field}"
                            if key in data: record[key] = data[key][idx[-1]:idx[-1]+1]
                        record[f"{name}_z"] = data[f"{name}_z"]
                    phase_records.append((Path(f"cycle_{int(cyc):04d}_phase_{int(ph)+1:03d}_recorded.npz"), record))
            window_ids = sorted(set(cycle_ids[cycle_ids >= report.get("window_start_cycle", int(cycle_ids.max()))]))
            window_ends = [np.flatnonzero(cycle_ids == k)[-1] for k in window_ids]
            fig, ax = plt.subplots(figsize=(8, 4))
            for name in names:
                pressure = np.asarray(data[f"{name}_P"])
                ax.plot(ht, trend(name, "P", np.mean(pressure, axis=-1)) / 1e5, label=name)
            ax.set(xlabel="Elapsed time [s]", ylabel="Mean bed pressure [bar]")
            ax.legend(title="Column")
            save(fig, "column_pressure_history")
            for name in names:
                prefix = safe_name(name)
                c = np.asarray(data[f"{name}_c"])
                q = np.asarray(data[f"{name}_q"])
                z = np.asarray(data[f"{name}_z"])
                y = c / np.maximum(c.sum(axis=1, keepdims=True), 1e-30)
                ce = np.stack([np.stack([trend(name, f"c_{end}_{j}", c[:, j, pos])
                                       for j in range(len(components))], axis=1)
                               for end, pos in (("left", 0), ("right", -1))], axis=-1)
                ye = ce / np.maximum(ce.sum(axis=1, keepdims=True), 1e-30)
                fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
                for index, component in enumerate(components):
                    axes[0].plot(ht, (ce if liquid else ye)[:, index, 0], label=component)
                    axes[1].plot(ht, (ce if liquid else ye)[:, index, -1], label=component)
                for ax, port in zip(axes, ("Left end", "Right end")):
                    ax.set(xlabel="Elapsed time [s]", title=f"{name}: {port}")
                    if not liquid:ax.set_ylim(-0.02,1.02)
                    ax.legend()
                axes[0].set_ylabel("Solute concentration [mol/m³ liquid]" if liquid else "Gas mole fraction [-]")
                save(fig, f"{prefix}_end_compositions")
                fig, axes = plt.subplots(1, 2, figsize=(10, 4))
                for index, component in enumerate(components):
                    for ax, values in ((axes[0], c), (axes[1], q)):
                        line, = ax.plot(z, values[selected_end, index], label=component)
                        ax.fill_between(z, values[window_ends, index].min(axis=0),
                                        values[window_ends, index].max(axis=0), color=line.get_color(), alpha=.15)
                axes[0].set(xlabel="Axial position [m]", ylabel="Solute concentration [mol/m³ liquid]" if liquid else "Gas concentration [mol/m³ gas]")
                axes[1].set(xlabel="Axial position [m]", ylabel="Total loading [mol/kg adsorbent]")
                for ax in axes:
                    ax.legend()
                fig.suptitle(f"Column {name}: selected end profiles; shading = available tail-cycle range")
                save(fig, f"{prefix}_final_profiles")
                fig, ax = plt.subplots(figsize=(8, 4))
                for field, label in (("Tg", "Liquid" if liquid else "Gas"), ("Ts", "Solid"), ("Tw", "Wall"), ("Tj", "Jacket/HTF")):
                    key = f"{name}_{field}"
                    if key in data:
                        ax.plot(ht, trend(name, field, np.mean(data[key], axis=-1)), label=label)
                ax.set(xlabel="Elapsed time [s]", ylabel="Mean temperature [K]", title=f"Column {name}")
                ax.legend()
                save(fig, f"{prefix}_temperature_history")
                if f"{name}_Tj" in data:
                    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
                    for field in ("Tg", "Ts", "Tw", "Tj"):
                        axes[0].plot(z, data[f"{name}_{field}"][selected_end], label=field)
                    active = trend(name, "htf_active", np.asarray(data[f"{name}_htf_active"])) > 0
                    for field, label in (("htf_inlet_temperature", "HTF inlet"), ("htf_outlet_temperature", "HTF outlet")):
                        axes[1].plot(ht, np.where(active, trend(name, field, data[f"{name}_{field}"]), np.nan), label=label)
                    axes[0].set(xlabel="Axial position [m]", ylabel="Temperature [K]", title="Selected cycle end")
                    axes[1].set(xlabel="Elapsed time [s]", ylabel="Temperature [K]", title="Active HTF stages")
                    for ax in axes: ax.legend()
                    save(fig, f"{prefix}_jacket_profiles")
                if selected_indices.size:
                    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
                    for cyc in (plot_cycles or [chosen]):
                        idx = np.flatnonzero(cycle_ids == cyc)
                        if not idx.size: continue
                        tc = t[idx] - cycle_starts.get(cyc, t[idx[0]])
                        for j, component in enumerate(components):
                            axes[0].plot(tc, (c if liquid else y)[idx, j, -1], label=f"{component}, cycle {cyc}")
                        for field in ("Tg", "Ts", "Tw", "Tj"):
                            if f"{name}_{field}" in data:
                                axes[1].plot(tc, np.mean(data[f"{name}_{field}"][idx], axis=-1), label=f"{field}, cycle {cyc}")
                    axes[0].set(xlabel="Time within cycle [s]", ylabel="Right-end concentration [mol/m³ liquid]" if liquid else "Right-end mole fraction [-]")
                    axes[1].set(xlabel="Time within cycle [s]", ylabel="Mean temperature [K]")
                    for ax in axes: ax.legend(fontsize=7)
                    save(fig, f"{prefix}_selected_cycle_history")
                if reference_values is not None:
                    pressure, temperature, composition, gas_constant = reference_values
                    reference = pressure * np.asarray(composition) / (gas_constant * temperature)
                    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
                    plotted = False
                    for index, component in enumerate(components):
                        if index >= reference.size or reference[index] <= 0:
                            continue
                        plotted = True
                        for ax, end in zip(axes, (0, -1)):
                            ax.plot(ht, ce[:, index, end] / reference[index], label=component)
                    if plotted:
                        for ax, port in zip(axes, ("Left", "Right")):
                            ax.set(xlabel="Elapsed time [s]", title=f"{name}: {port} end")
                            ax.legend()
                        axes[0].set_ylabel("C at bed end / fixed feed reference [-]")
                        fig.suptitle(f"Reference: {pressure / 1e5:g} bar, {temperature:g} K; "
                                     "not an instantaneous outlet ratio")
                        save(fig, f"{prefix}_normalised_end_concentrations")
                    else:
                        plt.close(fig)
    metric_path = output / "cycle_metrics.csv"
    if metric_path.is_file():
        with metric_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if rows:
            fig, ax = plt.subplots(figsize=(9, 5))
            cx = [int(row["cycle"]) for row in rows]
            errors = [float(row.get("css_state_error") or "nan") for row in rows]
            ax.semilogy(cx, np.maximum(errors, 1e-16), "o-", label="One-cycle state residual")
            periods = [json.loads(row.get("css_period_errors") or "{}") for row in rows]
            for lag in sorted({k for item in periods for k in item}):
                values = [item.get(lag, {}).get("css_state_error", np.nan) for item in periods]
                ax.semilogy(cx, np.maximum(np.asarray(values, dtype=float), 1e-16), "--", label=f"Lag {lag}")
            ax.axhline(1., color="black", linestyle=":", label="State criterion (other criteria also required)")
            if chosen is not None: ax.axvline(chosen, color="#aa4020", linestyle=":", label="Selected cycle")
            if report.get("window_start_cycle") is not None:
                ax.axvspan(report["window_start_cycle"], report["window_end_cycle"], color="grey", alpha=.1)
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            ax.set(xlabel="Completed cycle / full network schedule", ylabel="Tolerance-scaled CSS state error [-]")
            ax.legend(fontsize=8)
            save(fig, "css_convergence")
            groups = {"purity_recovery": [], "energy": []}
            for key in rows[0]:
                if any(term in key.lower() for term in ("purity", "recovery")):
                    groups["purity_recovery"].append(key)
                elif key.endswith("_J_cycle"):
                    groups["energy"].append(key)
            for stem, keys in groups.items():
                if not keys:
                    continue
                fig, ax = plt.subplots(figsize=(8, 4))
                for key in keys:
                    try:
                        values = [float(row.get(key) or "nan") for row in rows]
                    except (TypeError, ValueError):
                        continue
                    ax.plot([int(row.get("cycle", i + 1)) for i, row in enumerate(rows)], values,
                            marker="o", label=key)
                if ax.lines:
                    ax.set(xlabel="Completed cycle",
                           ylabel="Fraction [-]" if stem == "purity_recovery" else "Energy [J/cycle]")
                    ax.legend(fontsize=8)
                    save(fig, stem)
                else:
                    plt.close(fig)
    phase_files = sorted((output / "phases").glob("cycle_*_phase_*.npz"))
    prefixes = {f"cycle_{int(cyc):04d}" for cyc in plot_cycles}
    if phase_files:
        if not prefixes: prefixes = {phase_files[-1].name.split("_phase_", 1)[0]}
        phase_records = [(p, p) for p in phase_files if p.name.split("_phase_", 1)[0] in prefixes]
    if phase_records:
        phase_files = [item[0] for item in phase_records]
        if not prefixes: prefixes = {p.name.split("_phase_", 1)[0] for p in phase_files}
        cycle_prefix = ", ".join(sorted(prefixes))
        colours = plt.get_cmap("turbo")(np.linspace(0.05, .95, len(phase_files)))
        for name in names:
            fig, axes = plt.subplots(2, 2, figsize=(11, 8))
            for phase_index, ((phase_file, source), colour) in enumerate(zip(phase_records, colours), 1):
                with (np.load(source, allow_pickle=False) if isinstance(source, Path) else nullcontext(source)) as phase:
                    if f"{name}_c" not in phase:
                        continue
                    c = np.asarray(phase[f"{name}_c"])[-1]
                    q = np.asarray(phase[f"{name}_q"])[-1]
                    z = np.asarray(phase[f"{name}_z"])
                    y = c / np.maximum(c.sum(axis=0), 1e-30)
                    for component in range(len(components)):
                        style = ("-", "--", ":", "-.")[component % 4]
                        axes[0, 0].plot(z, (c if liquid else y)[component], color=colour, linestyle=style)
                        axes[0, 1].plot(z, q[component], color=colour, linestyle=style)
                    axes[1, 0].plot(z, phase[f"{name}_P"][-1] / 1e5, color=colour,
                                    label=phase_file.stem.split("_", 4)[1] + ":" + phase_file.stem.split("_", 4)[3])
                    axes[1, 1].plot(z, phase[f"{name}_Tg"][-1], color=colour)
                    axes[1, 1].plot(z, phase[f"{name}_Ts"][-1], color=colour, linestyle="--")
                    if f"{name}_Tj" in phase:
                        axes[1, 1].plot(z, phase[f"{name}_Tj"][-1], color=colour, linestyle=":")
            for ax in axes.flat:
                ax.set_xlabel("Axial position [m]")
            axes[0, 0].set_ylabel("Solute concentration [mol/m³ liquid]" if liquid else "Gas mole fraction [-]")
            axes[0, 1].set_ylabel("Total loading [mol/kg adsorbent]")
            axes[1, 0].set_ylabel("Pressure [bar]")
            axes[1, 1].set_ylabel("Temperature [K]")
            axes[1, 1].set_title("Gas: solid; solid phase: dashed; HTF: dotted", fontsize=9)
            axes[1, 0].legend(title="Phase", ncol=min(3, len(phase_files)), fontsize=7, loc="best")
            styles = "; ".join(f"{component}: {('solid', 'dashed', 'dotted', 'dash-dot')[i % 4]}"
                               for i, component in enumerate(components))
            fig.suptitle(f"Column {name}: end-of-phase profiles, {cycle_prefix.replace('_', ' ')}\n"
                         f"Component lines: {styles}")
            save(fig, f"{safe_name(name)}_phase_end_profiles")
    manifest = {
        "files": saved,
        "backend": "Matplotlib Agg",
        "generator": "tools/plot_results.py",
        "convergence_status": status,
        "selected_cycles": plot_cycles,
        "notes": [
            "Bed-end compositions are not necessarily flowing outlet compositions during closed or reversed-flow stages.",
            "Normalised curves use the fixed reference feed concentration; zero-feed components are omitted.",
            "Figures are drawn from the saved NPZ and CSV results; no equations are solved here.",
            "Compact whole-run histories use history_trends.csv; axial profiles use the retained diagnostic window.",
        ],
    }
    (plots / "plot_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return saved


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", nargs="?", default="results",
                        help="results directory written by the native run (default: results)")
    parser.add_argument("--formats", nargs="+", default=["png", "pdf"],
                        help="output formats (default: png pdf)")
    parser.add_argument("--dpi", type=int, default=160, help="raster resolution (default: 160)")
    parser.add_argument("--legacy-input", type=Path,
                        help="plot original column/component text outputs beside this input JSON")
    args = parser.parse_args(argv)
    if args.legacy_input:
        from plot_legacy_results import plot_legacy_results
        saved = plot_legacy_results(args.legacy_input, args.results, tuple(args.formats), args.dpi)
    else:
        saved = plot_results(args.results, tuple(args.formats), args.dpi)
    if not saved:
        raise SystemExit("No figures were produced")
    print(f"Wrote {len(saved)} plot file(s) to {Path(args.results) / 'plots'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
