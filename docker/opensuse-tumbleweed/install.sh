#!/bin/sh
# Build dependencies for the linux-opensuse-tumbleweed preset (default GCC).
set -eu

zypper --non-interactive install --no-recommends \
  git cmake ninja gcc gcc-c++ \
  suitesparse-devel lapack-devel python3-devel
zypper clean --all

python3 -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir nanobind numpy pandas plotly pytest
