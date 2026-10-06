# Build images

One directory per Linux distribution. `install.sh` is the complete list of packages needed to build Ruptura
there; the `Dockerfile` runs it on top of the official image, and the GitHub workflows run the same script inside
the same official image. Each image belongs to the CMake preset of the same distribution:

| Directory             | Base image                 | Preset                      | Compiler                           |
| --------------------- | -------------------------- | --------------------------- | ---------------------------------- |
| `ubuntu-24`           | `ubuntu:24.04`             | `linux-ubuntu-24`           | Clang 18, libc++                   |
| `ubuntu-26`           | `ubuntu:26.04`             | `linux-ubuntu-26`           | Clang 21, libc++                   |
| `debian-13`           | `debian:13`                | `linux-debian-13`           | Clang 19, libc++                   |
| `fedora-44`           | `fedora:44`                | `linux-fedora-44`           | GCC 16                             |
| `fedora-44`           | `fedora:44`                | `linux-clang`               | Clang 22, libc++                   |
| `rocky-9`             | `rockylinux/rockylinux:9`  | `linux-redhat-9`            | Clang 21, libc++ built from source |
| `rocky-10`            | `rockylinux/rockylinux:10` | `linux-redhat-10`           | Clang 21, libc++ built from source |
| `opensuse-leap-15.6`  | `opensuse/leap:15.6`       | `linux-opensuse-leap-15.6`  | Clang 19, libc++                   |
| `opensuse-tumbleweed` | `opensuse/tumbleweed`      | `linux-opensuse-tumbleweed` | GCC 16                             |
| `archlinux`           | `archlinux:base`           | `linux-archlinux`           | GCC 16 (x86_64 only)               |

Ruptura uses `<print>` and `<mdspan>`, so it needs GCC 16 or newer, or Clang with libc++ 17 or newer. Where the
distribution's default GCC is older, the image uses Clang with libc++ instead.

All images put CMake, Ninja, SuiteSparse (KLU), LAPACK, Python with its headers, and a virtual environment in
`/opt/venv` with nanobind, NumPy, pandas, Plotly and pytest. SUNDIALS and GoogleTest are downloaded by CMake.

## Running locally with Docker

From the repository root:

```
docker build -t ruptura-fedora-44 docker/fedora-44
docker run --rm -v "$PWD":/ruptura ruptura-fedora-44 cmake --workflow --preset linux-fedora-44
docker run --rm -v "$PWD":/ruptura ruptura-fedora-44 python -m pytest tests
```

The build goes to `build/<preset>` in your checkout, exactly as a native build would (on a Linux host the files are
owned by root unless you add `--user "$(id -u):$(id -g)"`). On Apple Silicon the images run natively as
`linux/arm64`; add `--platform linux/amd64` to `docker build` and `docker run` for x86_64.

## Running the GitHub workflows locally with act

[act](https://github.com/nektos/act) runs the Linux jobs of `.github/workflows` in Docker (0.2.89 works; 0.2.79 is
too old for the `node24` actions, `brew upgrade act`). The runner images are mapped in `.actrc`. For example:

```
act workflow_dispatch -W .github/workflows/test-matrix.yml -j linux --matrix preset:linux-fedora-44
act workflow_dispatch -W .github/workflows/test-matrix.yml -j linux    # every image
act workflow_dispatch -W .github/workflows/pull-request-checks.yml -j conda
act workflow_dispatch -W .github/workflows/conda-packages.yml -j build --matrix target-platform:linux-aarch64
```

act copies the working tree, including uncommitted and untracked files, so a job can pass locally while the same
job on GitHub fails on a file that was never committed. macOS jobs are skipped; `-P macos-15=-self-hosted` runs
them directly on a Mac instead. `--rm` in `.actrc` removes the containers of failed jobs; drop it to keep them
for inspection.

## Adding a distribution

1. Copy a directory, change `FROM` in the `Dockerfile` and the packages in `install.sh`.
2. Add a `linux-<distro>` configure, build, test and workflow preset to `CMakePresets.json`.
3. Add the image and preset to the `linux` matrix in `.github/workflows/test-matrix.yml`.
