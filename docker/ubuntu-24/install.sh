#!/bin/sh
# Build dependencies for the linux-ubuntu-24 preset (Clang 18 with libc++; GCC 13 lacks <print> and <mdspan>).
set -eu

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
  ca-certificates git cmake ninja-build clang libc++-dev libc++abi-dev \
  libsuitesparse-dev liblapack-dev python3-dev python3-venv
rm -rf /var/lib/apt/lists/*

python3 -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir nanobind numpy pandas plotly pytest
