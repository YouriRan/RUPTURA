#!/bin/sh
# Build dependencies for the linux-ubuntu-26 preset: Clang 21 with libc++ (the default GCC 15 has no <mdspan>;
# g++-16 is only an experimental snapshot).
set -eu

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
  ca-certificates git cmake ninja-build clang libc++-dev libc++abi-dev \
  libsuitesparse-dev liblapack-dev python3-dev python3-venv
rm -rf /var/lib/apt/lists/*

python3 -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir nanobind numpy pandas plotly pytest
