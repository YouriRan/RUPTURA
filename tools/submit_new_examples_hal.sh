#!/usr/bin/env bash
# Submit the selected examples as one array; no concurrency cap by default.
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec python3 "$root/tools/hal_examples.py" submit "$@"
