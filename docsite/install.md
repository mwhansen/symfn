# Installing

symfn is a compiled extension module with a pure-Python layer on top. It has no
Python dependencies at all, and it does not require, import, or know about Sage.

```console
$ pip install symfn
```

The wheel is built `abi3`, so one artifact per platform serves CPython 3.9 and
later.

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

## Checking the install

```pycon
>>> import symfn
>>> symfn.__version__
'0.1.0'
>>> symfn.s([2, 1]) * symfn.s([1])
s[2,1,1] + s[2,2] + s[3,1]
```

## Using it from Sage

Sage is a *consumer* of this library, on the far side of the boundary — the
wheel itself contains no Sage code and gains nothing from Sage being present.
The adapter that lets Sage's own symmetric-function classes compute through
symfn lives in the repository under `scripts/`, not in the wheel; installing
symfn does not change Sage's behavior on its own.
