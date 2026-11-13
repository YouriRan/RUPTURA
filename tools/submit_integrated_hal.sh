#!/usr/bin/env bash
# Compatibility entry point for the original two campaign names.
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
selection="${1:-cycles}"
if (($#)); then shift; fi
case "$selection" in
  cycles) selection=mfi_cycles ;;
  breakthrough) selection=mfi_breakthrough ;;
esac
exec "$root/tools/submit_new_examples_hal.sh" "$selection" "$@"
