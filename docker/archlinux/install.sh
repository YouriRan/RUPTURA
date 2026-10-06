#!/bin/sh
# Build dependencies for the linux-archlinux preset (default GCC). The archlinux image is x86_64 only.
set -eu

pacman -Syu --noconfirm --needed --disable-sandbox git cmake ninja gcc suitesparse lapack python
pacman -Scc --noconfirm

python -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir nanobind numpy pandas plotly pytest
