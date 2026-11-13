#!/bin/sh
# Build dependencies for the linux-opensuse-leap-15.6 preset: Clang 19 with libc++ (gcc15 has no <mdspan>; gcc is
# still needed for the C runtime start files). Python 3.11.
set -eu

zypper --non-interactive install --no-recommends \
  git cmake ninja gcc clang19 libc++-devel libc++abi-devel \
  suitesparse-devel lapack-devel python311-devel
zypper clean --all

python3.11 -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir nanobind numpy pandas plotly pytest
