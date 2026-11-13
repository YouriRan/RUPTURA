#!/usr/bin/env bash
#SBATCH --partition=serial
#SBATCH --exclude=c171
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname -- "${BASH_SOURCE[0]}")}"
exec ./run "$@"
