# Installing

symfn is a compiled extension module with a pure-Python layer on top. It has no
Python dependencies at all.

```console
$ pip install symfn
```

The wheel is built `abi3`, so one artifact per platform serves CPython 3.9 and
later.

## Platforms

A prebuilt wheel covers Linux (glibc and musl) on x86\_64, aarch64, i686,
armv7l, ppc64le and s390x; macOS on x86\_64 and arm64; and Windows on amd64,
win32 and arm64.

Anywhere else (FreeBSD, illumos, riscv64, and any other target Rust
supports) needs to install from the source distribution and as such
needs a Rust toolchain installed. That build needs no network: the sdist carries
its four Rust dependencies vendored inside it.

## Building from source

Building needs a Rust toolchain (1.87 or later) and
[maturin](https://www.maturin.rs/):

```console
$ pip install maturin
$ maturin build --release
```

or, to build and install into the active environment in one step:

```console
$ maturin develop --release
```

The `--release` matters for more than speed: the release profile keeps overflow
checks on, which is part of the correctness surface rather than a debug aid.
