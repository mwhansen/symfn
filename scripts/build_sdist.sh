#!/bin/sh
# Build the source distribution, with the wheel build's dependencies vendored
# into it so it compiles with the network off.
#
#     scripts/build_sdist.sh [outdir]        # default outdir: dist/
#
# Why vendoring is part of building the sdist rather than a separate courtesy:
# the default crate has zero dependencies and builds offline already, but the
# `python` feature pulls PyO3, num-bigint, num-rational and num-traits, and an
# offline source build is exactly the configuration distro packagers use
# (docs/support-tiers.md, Tier 2). Sage itself never builds Rust — it consumes
# the prebuilt wheels (docs/sage-packaging-audit.md) — so this artifact is for
# Debian, conda-forge, Gentoo and nix, who build from source on principle.
#
# Three things make the result offline: `vendor/` holds the crate sources,
# `.cargo/config.toml` redirects crates-io at them, and `Cargo.lock` pins what
# was vendored. The lockfile is tracked, so `--locked` is what makes the
# vendored set and the published lockfile the same set rather than whatever
# resolved today.
#
# **`.cargo/config.toml` is tracked and carries the macOS linker flags a bare
# `cargo build --features python` needs**, so the redirect is appended to it
# and the original is put back on the way out. An earlier draft wrote the file
# outright, which deleted those flags and broke every local PyO3 build until
# the next `git checkout`. The trap on EXIT is what makes an interrupted run
# restore it too.
#
# `scripts/check_sdist_offline.sh` is the assertion that this worked.

set -e

here=$(dirname "$0")
root=$(cd "$here/.." && pwd)
out=${1:-$root/dist}

cd "$root"

restore() {
	rm -rf "$root/vendor"
	if [ -f "$root/.cargo/config.toml.sdist-backup" ]; then
		mv "$root/.cargo/config.toml.sdist-backup" "$root/.cargo/config.toml"
	fi
}
trap restore EXIT

# A leftover vendor directory from an interrupted run would be vendored *into*
# the next one.
rm -rf vendor

printf '== sdist: vendoring the python feature dependencies\n'
cargo vendor --locked --versioned-dirs vendor >/dev/null

mkdir -p .cargo
cp .cargo/config.toml .cargo/config.toml.sdist-backup
cat >>.cargo/config.toml <<'EOF'

# Appended by scripts/build_sdist.sh; present only inside the sdist.
#
# `directory` is relative to the project root -- the parent of this file's
# directory -- so the sdist is relocatable: unpack it anywhere and cargo reads
# the crate sources beside it rather than reaching for crates.io.
[source.crates-io]
replace-with = "vendored-sources"

[source.vendored-sources]
directory = "vendor"
EOF

printf '== sdist: maturin sdist\n'
maturin sdist --out "$out"

printf '== sdist: built in %s\n' "$out"
ls -l "$out"
