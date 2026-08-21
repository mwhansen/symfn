#!/bin/sh
# Build both references and bundle them into one downloadable archive.
#
#     scripts/build_docs.sh [outdir]        # default outdir: dist/
#
# Two references, because there are two audiences and neither one's manual
# covers the other's questions. `cargo doc` is the Rust API — the traits, the
# generic parameters, the module docs that carry the convention sections. The
# Sphinx site is the Python surface, generated from the objects in the built
# wheel so `help()` and the website cannot disagree. A landing page at the root
# of the bundle points at both.
#
# Why a bundle at all: Read the Docs publishes the Python half and nothing
# publishes the Rust half, so until both have a public home a tester who
# downloaded a wheel has no reference to read beside it. This archive is
# attached to every GitHub Release for exactly that stretch, and it is
# self-contained — open `index.html` from a local filesystem, no server.
#
# `cargo doc` runs with warnings denied. CLAUDE.md states the tree keeps
# `cargo doc --no-deps --all-features` silent; this is what makes that a claim
# that can be false.

set -e

here=$(dirname "$0")
root=$(cd "$here/.." && pwd)
out=${1:-$root/dist}

cd "$root"

version=$(python3 -c "
import re, pathlib
text = pathlib.Path('Cargo.toml').read_text()
print(re.search(r'^version = \"([^\"]+)\"', text, re.M).group(1))
")

printf '== docs: cargo doc --no-deps --all-features\n'
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features

printf '== docs: staging the extension module\n'
# The Sphinx half documents the objects in the built module rather than a
# parsed copy of the source, which is what keeps `help()` and the website from
# disagreeing -- so the extension has to exist and sit inside the package
# before the site can be generated. `docsite/conf.py` only puts the source tree
# on `sys.path` when it finds one there.
cargo build --features python
python3 -c "
import sys
sys.path.insert(0, 'scripts')
from check_convenience_docs import stage
if stage() is None:
    raise SystemExit('no extension module after cargo build --features python')
"

printf '== docs: sphinx\n'
# `-W` for the same reason `cargo doc` denies warnings: the house docstring
# form is translated to RST by a hook in docsite/conf.py that a new
# construction can escape, and a warning is how that gets noticed.
python3 -m sphinx -b html -W docsite docsite/_build/html >/dev/null

stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
bundle=$stage/symfn-docs-$version
mkdir -p "$bundle"

cp -R docsite/_build/html "$bundle/python"
cp -R target/doc "$bundle/rust"

# `cargo doc` writes its landing page one level down, under the crate name, and
# the crate list at the root is a redirect that assumes it was served from the
# same path. A reader opening the archive should not have to know that.
cat >"$bundle/index.html" <<EOF
<!doctype html>
<meta charset="utf-8">
<title>symfn $version — references</title>
<style>
  body { font: 16px/1.6 system-ui, sans-serif; margin: 4rem auto; max-width: 34rem; padding: 0 1rem; }
  h1 { font-size: 1.5rem; }
  li { margin: .75rem 0; }
</style>
<h1>symfn $version</h1>
<p>Two references, one for each surface.</p>
<ul>
  <li><a href="python/index.html">Python</a> — the convenience layer, the
      families, and every contract entry point.</li>
  <li><a href="rust/symfn/index.html">Rust</a> — the crate API.</li>
</ul>
<p>Both are generated from the code they describe. This bundle is
   self-contained: no network, no server.</p>
EOF

mkdir -p "$out"
tar czf "$out/symfn-docs-$version.tar.gz" -C "$stage" "symfn-docs-$version"

printf '== docs: bundled in %s\n' "$out"
ls -l "$out/symfn-docs-$version.tar.gz"
