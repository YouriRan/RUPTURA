#!/bin/sh
# Build dependencies for the linux-redhat-9 preset: Clang from AppStream with libc++ built from the matching LLVM
# sources, because gcc-toolset-15 has no <mdspan> (switch to gcc-toolset-16 once RHEL ships it). CMake and Ninja
# come from pip, since the system CMake is older than 3.28 and the system Python older than 3.11.
set -eu

dnf install -y 'dnf-command(config-manager)'
dnf config-manager --set-enabled crb
dnf install -y epel-release
dnf install -y --setopt=install_weak_deps=False \
  git xz gcc clang \
  suitesparse-devel lapack-devel python3.12 python3.12-devel
dnf clean all

python3.12 -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir cmake ninja nanobind numpy pandas plotly pytest

version=$(clang -dumpversion)
curl -sSL "https://github.com/llvm/llvm-project/releases/download/llvmorg-${version}/llvm-project-${version}.src.tar.xz" \
  | tar -xJ -C /tmp
/opt/venv/bin/cmake -G Ninja -S "/tmp/llvm-project-${version}.src/runtimes" -B /tmp/libcxx \
  -DCMAKE_MAKE_PROGRAM=/opt/venv/bin/ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_INSTALL_PREFIX=/usr \
  -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind" -DLLVM_ENABLE_PER_TARGET_RUNTIME_DIR=OFF -DLLVM_LIBDIR_SUFFIX=64 \
  -DLIBCXX_INCLUDE_TESTS=OFF -DLIBCXX_INCLUDE_BENCHMARKS=OFF -DLIBCXXABI_INCLUDE_TESTS=OFF -DLIBUNWIND_INCLUDE_TESTS=OFF
/opt/venv/bin/ninja -C /tmp/libcxx install-cxx install-cxxabi install-unwind
rm -rf "/tmp/llvm-project-${version}.src" /tmp/libcxx
ldconfig
