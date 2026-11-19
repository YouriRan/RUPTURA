"""TikZ drawing of the columns, buffers, valves, and compressors of a swing-adsorption input with `Columns`.

    python -m ruptura.flowsheet simulation.json [phase] > flowsheet.tex

draws every connection of the cycle; with a phase (its name, or its number from 1) the connections of that phase are
black with the direction of From to To, the others gray, and the ports of the phase are at the column ends. The
picture needs \\usepackage{tikz} only. Columns stand upright with the feed end (z = 0) at the bottom.
"""

from __future__ import annotations

import json
import sys

from .webapp.config import canonicalize

WIDTH, HEIGHT, PITCH, LANE = 0.8, 3.0, 2.6, 0.7  # of a column, between columns, and between pipe runs, in cm

SYMBOLS = r"""  % Symbols: a valve, a check valve (filled towards its outlet), and a compressor, each with its outlet at +x.
  \tikzset{
    valve/.pic={\draw[fill=white] (-0.16,-0.1) -- (0.16,0.1) -- (0.16,-0.1) -- (-0.16,0.1) -- cycle;},
    check valve/.pic={\draw[fill=white] (-0.16,-0.1) -- (0,0) -- (-0.16,0.1) -- cycle;
                      \draw[fill] (0.16,-0.1) -- (0,0) -- (0.16,0.1) -- cycle;},
    compressor/.pic={\draw[fill=white] (0,0) circle (0.18); \draw (-0.09,-0.156) -- (0.18,0) -- (-0.09,0.156);}}
"""


def _text(value):
    """A name as LaTeX text."""
    return "".join("\\" + character if character in "_%&#$" else character for character in str(value))


def flowsheet_tikz(data, phase=None):
    """The tikzpicture of the input data (a dict as read from simulation.json); see the module description."""
    data = canonicalize(data)
    columns = list(data["Columns"])
    buffers = list(data.get("Buffers", {}))
    phases = data["SwingAdsorptionPhases"]
    if isinstance(phase, str):
        phase = [item.get("Name", "").lower() for item in phases].index(phase.lower()) + 1
    shown = phases[phase - 1] if phase else {}

    # Every connection of the cycle once, by its ends and its name; those of the shown phase keep their direction. An
    # end is (column or buffer index, "FeedEnd", "ProductEnd", or "Buffer").
    connections = {}
    for item in phases:
        for connection in item.get("Connections", []):
            ends = []
            for name in (connection["From"], connection["To"]):
                owner, _, end = name.partition(".")
                kind = {"feedend": "FeedEnd", "productend": "ProductEnd", "": "Buffer"}[end.lower()]
                names = buffers if kind == "Buffer" else columns
                ends.append(([other.lower() for other in names].index(owner.lower()), kind))
            ends = tuple(ends)
            key = (frozenset(ends), connection.get("Name", ""))
            if key not in connections or item is shown:
                connections[key] = {**connection, "ends": ends, "active": item is shown}

    # Pipe runs: those that touch only product ends and buffers lie above the columns, those of feed ends below, and a
    # connection from a feed end to anything above takes a run on both sides and a riser to the right of the columns.
    above = [value for value in connections.values() if any(end[1] != "FeedEnd" for end in value["ends"])]
    above.sort(key=lambda value: any(end[1] == "Buffer" for end in value["ends"]))  # the runs to buffers on top
    below = [value for value in connections.values() if any(end[1] == "FeedEnd" for end in value["ends"])]
    gap = 1.3 if phase else 0.6  # between a column end and the first run: room for the port of a shown phase
    buffer_y = HEIGHT + gap + LANE * len(above) + 0.5
    right = PITCH * max(len(columns) - 1, len(buffers) - 0.5) + 1.2
    position = lambda end: (PITCH * end[0], {"FeedEnd": 0.0, "ProductEnd": HEIGHT}[end[1]]) if end[1] != "Buffer" else (
        PITCH * (end[0] + 0.5), buffer_y)

    lines = [r"\begin{tikzpicture}[line cap=round, line join=round, font=\small]", SYMBOLS.rstrip("\n")]
    for index, name in enumerate(columns):
        x = PITCH * index
        lines.append(rf"  \draw[thick] ({x - WIDTH / 2:.2f},0) rectangle ({x + WIDTH / 2:.2f},{HEIGHT:.2f});"
                     rf" \node at ({x:.2f},{HEIGHT / 2:.2f}) {{{_text(name)}}};")
    for index, name in enumerate(buffers):
        x, y = position((index, "Buffer"))
        lines.append(rf"  \node[draw, thick, rounded corners=0.2cm, minimum width=1.2cm, minimum height=0.8cm,"
                     rf" anchor=south] at ({x:.2f},{y:.2f}) {{{_text(name)}}};")

    risers = 0
    for value in connections.values():
        (x0, y0), (x1, y1) = position(value["ends"][0]), position(value["ends"][1])
        kinds = [end[1] for end in value["ends"]]
        high = HEIGHT + gap + LANE * above.index(value) if value in above else None
        low = -gap - LANE * below.index(value) if value in below else None
        run = low if high is None else high  # the run that carries the symbol
        if high is not None and low is not None:
            riser = right + LANE * risers
            risers += 1
            turn = lambda kind: low if kind == "FeedEnd" else high
            path = [(x0, y0), (x0, turn(kinds[0])), (riser, turn(kinds[0])), (riser, turn(kinds[1])),
                    (x1, turn(kinds[1])), (x1, y1)]
            run, start, end = turn(kinds[0]), x0, riser
        else:
            path, start, end = [(x0, y0), (x0, run), (x1, run), (x1, y1)], x0, x1
        symbol = ("compressor" if value.get("FlowModel", "").lower() == "compressor"
                  else "check valve" if value.get("CheckValve") else "valve")
        style = "" if value["active"] or not phase else "[gray!60]"
        arrow = "[->]" if value["active"] else style
        points = " -- ".join(f"({x:.2f},{y:.2f})" for x, y in path)
        lines.append(rf"  \draw{arrow or style} {points};")
        lines.append(rf"  \pic{style}[xscale={1 if end >= start else -1}] at ({(start + end) / 2:.2f},{run:.2f})"
                     rf" {{{symbol}}}; \node{style}[above=0.08cm] at ({(start + end) / 2:.2f},{run:.2f})"
                     rf" {{\footnotesize {_text(value.get('Name', ''))}}};")

    # The ports of the shown phase: its own FeedEnd and ProductEnd, or those of its Columns entry for a column.
    for index, name in enumerate(columns if phase else []):
        own = next((value for key, value in shown.get("Columns", {}).items() if key.lower() == name.lower()), {})
        for key, y, direction in (("FeedEnd", 0.0, -1), ("ProductEnd", HEIGHT, 1)):
            port = own.get(key, shown.get(key))
            if port is None or port.get("Mode", "").lower() == "closed":
                continue
            # A pressure port passes gas either way; a velocity port carries its velocity or molar flux inwards.
            x, pressure = PITCH * index, port["Mode"].lower() == "pressure"
            label = (f"{port['Pressure'] / 1e5:g} bar" if pressure else f"{port['Velocity']:g} m/s"
                     if "Velocity" in port else f"{port['MolarFlux']:g} mol/(m$^2$ s)")
            tip = "<->" if pressure else "<-"
            lines.append(rf"  \draw[{tip}, thick] ({x:.2f},{y:.2f}) -- ({x:.2f},{y + 0.45 * direction:.2f})"
                         rf" node[{'above' if direction > 0 else 'below'}] {{\footnotesize {label}}};")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8") as file:
        argument = sys.argv[2] if len(sys.argv) > 2 else None
        sys.stdout.write(flowsheet_tikz(json.load(file), int(argument) if argument and argument.isdigit() else argument))
