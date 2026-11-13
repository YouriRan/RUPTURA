#!/usr/bin/env bash
# SLURM may copy this script to its spool directory. Resolve code from the manifest.
set -euo pipefail
manifest="${1:?A frozen run manifest is required}"
index="${SLURM_ARRAY_TASK_ID:?Run this through submit_new_examples_hal.sh}"
exec python3 - "$manifest" "$index" <<'PY'
import json
import os
from pathlib import Path
import sys

manifest = Path(sys.argv[1]).resolve()
project = Path(json.loads(manifest.read_text())["project_root"])
os.execv(sys.executable, [sys.executable, str(project / "tools/hal_examples.py"),
                         "run", str(manifest), "--index", sys.argv[2]])
PY
