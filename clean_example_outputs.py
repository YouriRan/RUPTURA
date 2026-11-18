#!/usr/bin/env python3
"""Remove Ruptura output files from the examples tree."""

from __future__ import annotations

import argparse
from pathlib import Path


OUTPUT_NAMES = {"column.data", "failed_state.json", "restart_data.bin", "run", "column.json", "breakthrough.ipynb", "nohup.out", 'newplot.png', "mixture.ipynb"}
OUTPUT_DIR_NAMES = {"output", "density_grids"}


def output_candidates(examples: Path) -> set[Path]:
    candidates = {
        path
        for path in examples.rglob("*")
        if path.is_file()
        and (path.name in OUTPUT_NAMES or path.name.startswith("component_") and path.suffix == ".data")
    }
    for directory in examples.rglob("*"):
        if directory.is_dir() and directory.name in OUTPUT_DIR_NAMES:
            candidates.update(path for path in directory.rglob("*") if path.is_file())
    return candidates


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print files that would be removed without changing anything",
    )
    parser.add_argument(
        "--examples-dir",
        type=Path,
        help="examples directory (defaults to examples/ beside this script)",
    )
    args = parser.parse_args()

    repository = Path(__file__).resolve().parent
    examples = (args.examples_dir or repository / "examples").resolve()
    if not examples.is_dir():
        parser.error(f"examples directory does not exist: {examples}")

    candidates = sorted(output_candidates(examples))

    action = "Would remove" if args.dry_run else "Removed"
    for path in candidates:
        print(path.relative_to(repository) if path.is_relative_to(repository) else path)
        if not args.dry_run:
            path.unlink()

    if not args.dry_run:
        for directory in sorted(
            (path for path in examples.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts),
            reverse=True,
        ):
            if directory.name in OUTPUT_DIR_NAMES:
                try:
                    directory.rmdir()
                except OSError:
                    pass

    print(f"{action} {len(candidates)} file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
