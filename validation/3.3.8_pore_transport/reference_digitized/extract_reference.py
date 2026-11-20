#!/usr/bin/env python3
"""Write the Goeppert et al. (2014) breakthrough data used by Auton et al. (2024) to one CSV file per particle size.

The input is `Figure_8_9_data.m` from `Intra-Particle-codes.zip` in https://github.com/aguareles/Intra-particle-diffusion,
which holds (time in h, outlet CO2 in ppm) as digitized by A. Valverde from Goeppert et al. The file is only parsed as
text, never executed.

    python extract_reference.py path/to/Figure_8_9_data.m
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Data function in the MATLAB file -> output file, with the sieve fraction and the radius Auton et al. assign to it.
FRACTIONS = {
    "fun_data_0o25": "goeppert2014_below_0p25mm_R1p25e-4m.csv",
    "fun_data_0o25_0o5": "goeppert2014_0p25_to_0p5mm_R1p875e-4m.csv",
    "fun_data_0o5_1o7": "goeppert2014_0p5_to_1p7mm_R5p5e-4m.csv",
    "fun_data_1o7": "goeppert2014_above_1p7mm_R8p5e-4m.csv",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("matlab_file", type=Path)
    arguments = parser.parse_args()

    text = arguments.matlab_file.read_text(errors="replace")
    for function_name, file_name in FRACTIONS.items():
        block = re.search(rf"function\s+\[?\s*Data\s*\]?\s*=\s*{function_name}\b.*?Data\s*=\s*\[(.*?)\]", text, re.S)
        if block is None:
            raise SystemExit(f"{function_name} not found in {arguments.matlab_file}")
        rows = [line.split() for line in block.group(1).splitlines() if line.strip()]
        (ROOT / file_name).write_text("time_h,outlet_CO2_ppm\n" + "".join(f"{t},{ppm}\n" for t, ppm in rows))
        print(f"{file_name}: {len(rows)} points")


if __name__ == "__main__":
    main()
