#!/usr/bin/env python3
"""Pack the repository for a full rerun elsewhere (validation/README.md, "Full rerun on a cluster node"), without
any RUPTURA output.

    python validation/pack.py                      # writes ruptura-rerun-<date>.tar.gz next to the repository
    python validation/pack.py -o /path/to/out.tar.gz
    python validation/pack.py --list               # print what would be packed, and stop

The working tree is packed, not a git checkout, because the paper is a separate repository (git-ignored here): the
sources, examples, every input, generator, runner, and notebook, and the paper, whose figure folders the notebooks
write to. Left out:

  * RUPTURA output, by file name, wherever it lies: component_*.data, column.data, column.json, phases.data,
    cycles.data, buffers.data, connections.data, failed_state.json, fitted.json and *_fitted.json (multicolumn runs
    write them into one subdirectory per column), and the run logs run.log, ruptura.log, run.sh.log, nohup.out, and
    rerun_results.csv. Other .data files are inputs and are kept, such as the MPD distributions
    validation/3.2.2_mpd/sharma/mpd/particle_distribution_*.data. The restart files of the Dantas runs (precondition/column.json) are
    written again by their run.sh.
  * .git (also that of the paper), the build trees (build* at the top), .worktrees, invruptura, tmp, simulations (the
    Ruptura Lab records), __pycache__, .ipynb_checkpoints, *.egg-info, compiled modules (*.so, *.dylib), .DS_Store.
  * validation/3.2.2_mpd/sharma/tmmc: the RASPA3 TMMC data behind the MPD distributions (11 GB, not a RUPTURA job).

Symlinks, such as the old case paths at the repository root, are packed as links. Pack once no run is writing.
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import tarfile
from datetime import date
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
OUTPUT_PATTERNS = (
    "component_*.data", "column.data", "column.json", "phases.data", "cycles.data", "buffers.data", "connections.data",
    "failed_state.json", "fitted.json", "*_fitted.json",
    "run.log", "ruptura.log", "run.sh.log", "nohup.out", "rerun_results.csv",
)
JUNK_PATTERNS = ("*.so", "*.dylib", ".DS_Store")
EXCLUDED_ANYWHERE = (".git", "__pycache__", ".ipynb_checkpoints")
EXCLUDED_AT_TOP = (".worktrees", "invruptura", "tmp", "simulations")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", type=Path,
                        default=REPOSITORY.parent / f"{REPOSITORY.name}-rerun-{date.today():%Y%m%d}.tar.gz")
    parser.add_argument("--list", action="store_true", help="print the packed paths and stop")
    arguments = parser.parse_args()

    # Walk the tree without following symlinks; a symlink, to a file or a directory, is packed as the link itself.
    members = []
    for current, directories, files in os.walk(REPOSITORY):
        relative = Path(current).relative_to(REPOSITORY)
        kept = []
        for name in sorted(directories):
            path = Path(current) / name
            if (name in EXCLUDED_ANYWHERE or name.endswith(".egg-info")
                    or (relative == Path(".") and (name.startswith("build") or name in EXCLUDED_AT_TOP))
                    or (relative / name).as_posix().endswith("sharma/tmmc")):
                continue
            if path.is_symlink():
                members.append(path)
            else:
                kept.append(name)
        directories[:] = kept
        members.extend(Path(current) / name for name in sorted(files)
                       if not any(fnmatch.fnmatch(name, pattern) for pattern in OUTPUT_PATTERNS + JUNK_PATTERNS))

    size = sum(path.lstat().st_size for path in members)
    inputs = sum(path.name == "simulation.json" for path in members)
    if arguments.list:
        for path in members:
            print(path.relative_to(REPOSITORY.parent))
        print(f"{len(members)} files and links, {size / 1e6:.0f} MB before compression, {inputs} simulation.json")
        return

    with tarfile.open(arguments.output, "w:gz", compresslevel=6) as archive:
        for path in members:
            archive.add(path, arcname=path.relative_to(REPOSITORY.parent).as_posix(), recursive=False)
    print(f"{arguments.output}: {len(members)} files and links, {size / 1e6:.0f} MB before compression, "
          f"{arguments.output.stat().st_size / 1e6:.0f} MB compressed, {inputs} simulation.json")


if __name__ == "__main__":
    main()
