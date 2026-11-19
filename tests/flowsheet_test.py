"""The TikZ flowsheet of an input with several columns."""

import json
from pathlib import Path

from ruptura.flowsheet import flowsheet_tikz

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "PSA-H2-activated-carbon-two-beds" / "simulation.json"


def test_flowsheet_of_the_two_bed_example():
    data = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    picture = flowsheet_tikz(data)
    assert picture.startswith(r"\begin{tikzpicture}") and picture.endswith("\\end{tikzpicture}\n")
    # Two columns and the tank; the product valve of each bed (check valves), the purge, and the equalization, each
    # once although both halves of the cycle use them.
    assert all(f"{{{name}}};" in picture for name in ("A", "B", "Tank"))
    assert picture.count("{check valve};") == 2 and picture.count("{valve};") == 2
    assert "gray" not in picture and "->" not in picture and " bar" not in picture

    # The adsorption of A: its product and purge valves carry gas, the others are gray, and the feed ends hold a
    # velocity port and a pressure port.
    shown = flowsheet_tikz(data, "adsorption a, purge b")
    assert shown == flowsheet_tikz(data, 2)
    assert shown.count(r"\draw[->]") == 2 and shown.count(r"\draw[gray!60]") == 2
    assert r"\draw[<-, thick]" in shown and "0.03 m/s" in shown
    assert r"\draw[<->, thick]" in shown and "1 bar" in shown


def test_flowsheet_routes_feed_ends_buffers_and_compressors_and_escapes_names():
    data = {"columns": ["A", "B"], "buffers": {"Tank_1": {}}, "swingadsorptionphases": [
        {"Name": "Recycle", "connections": [
            {"Name": "K1", "from": "a.feedend", "to": "tank_1", "FlowModel": "Compressor"},
            {"Name": "V1", "From": "A.FeedEnd", "To": "B.FeedEnd"}]}]}
    picture = flowsheet_tikz(data, 1)
    assert r"Tank\_1" in picture and "{compressor};" in picture
    # The compressor line leaves the feed end downwards, rises beside the columns, and ends at the tank: six points.
    compressor = next(line for line in picture.splitlines() if line.startswith(r"  \draw[->]"))
    assert compressor.count("--") == 5
