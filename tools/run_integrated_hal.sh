#!/usr/bin/env bash
set -euo pipefail
manifest="${1:?Case manifest is required}"
engine="${2:?Executable path is required}"
shift 2
mapfile -t inputs < "$manifest"
input="${inputs[${SLURM_ARRAY_TASK_ID:?Run this through submit_integrated_hal.sh}]}"
cd -- "$(dirname -- "$input")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
printf 'Case %s: %s\n' "$SLURM_ARRAY_TASK_ID" "$input"
exec srun "$engine" "$input" "$@"
