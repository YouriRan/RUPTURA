#!/bin/sh
# Build dependencies for the linux-fedora-44 preset (default GCC) and the linux-clang preset (Clang with libc++).
set -eu

dnf install -y --setopt=install_weak_deps=False \
  git cmake ninja-build gcc gcc-c++ clang libcxx-devel libcxxabi-devel \
  suitesparse-devel lapack-devel python3-devel
dnf clean all

python3 -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir nanobind numpy pandas plotly pytest
