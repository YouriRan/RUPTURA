#!/bin/sh
# Build dependencies for the linux-redhat-10 preset: Clang from AppStream with libc++ built from the matching LLVM
# sources, because GCC 14 and gcc-toolset-15 have no <mdspan> (switch to gcc-toolset-16 once RHEL ships it).
set -eu

dnf install -y 'dnf-command(config-manager)'
dnf config-manager --set-enabled crb
dnf install -y epel-release
dnf install -y --setopt=install_weak_deps=False \
  git xz cmake ninja-build gcc clang \
  suitesparse-devel lapack-devel python3-devel
dnf clean all

python3 -m venv /opt/venv
/opt/venv/bin/pip install --no-cache-dir nanobind numpy pandas plotly pytest

version=$(clang -dumpversion)
curl -sSL "https://github.com/llvm/llvm-project/releases/download/llvmorg-${version}/llvm-project-${version}.src.tar.xz" \
  | tar -xJ -C /tmp
cmake -G Ninja -S "/tmp/llvm-project-${version}.src/runtimes" -B /tmp/libcxx -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_INSTALL_PREFIX=/usr \
  -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind" -DLLVM_ENABLE_PER_TARGET_RUNTIME_DIR=OFF -DLLVM_LIBDIR_SUFFIX=64 \
  -DLIBCXX_INCLUDE_TESTS=OFF -DLIBCXX_INCLUDE_BENCHMARKS=OFF -DLIBCXXABI_INCLUDE_TESTS=OFF -DLIBUNWIND_INCLUDE_TESTS=OFF
ninja -C /tmp/libcxx install-cxx install-cxxabi install-unwind
rm -rf "/tmp/llvm-project-${version}.src" /tmp/libcxx
ldconfig
