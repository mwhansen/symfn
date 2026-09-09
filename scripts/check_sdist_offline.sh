#!/bin/sh
# Fail unless the source distribution builds with the network off.
#
#     scripts/build_sdist.sh
#     scripts/check_sdist_offline.sh [sdist]   # default: dist/symfn-*.tar.gz
#
# The proposition: **a distro packager who unpacks the sdist on a build machine
# with no network gets a working wheel** (README.md, Tier 2). That
# is the configuration Debian, conda-forge, Gentoo and nix build in, and it is
# the one nothing else in the tree exercises -- every other Python job installs
# from the working tree with crates.io reachable.
#
# Offline is asserted rather than simulated. Cutting the runner's network would
# test the runner; `CARGO_NET_OFFLINE=true` makes cargo *refuse* to fetch, so a
# dependency the vendoring missed is a hard error naming the crate rather than
# a silent download. `--no-index` does the same for pip, and `--no-build-isolation`
# keeps pip from reaching for maturin -- which is why maturin has to be
# installed already.

set -e

here=$(dirname "$0")
root=$(cd "$here/.." && pwd)

sdist=$1
if [ -z "$sdist" ]; then
	sdist=$(ls "$root"/dist/symfn-*.tar.gz 2>/dev/null | head -1)
fi
if [ -z "$sdist" ] || [ ! -f "$sdist" ]; then
	echo "no sdist; build one: scripts/build_sdist.sh"
	exit 1
fi

if ! python3 -c "import maturin" 2>/dev/null && ! command -v maturin >/dev/null 2>&1; then
	echo "maturin is not installed, and --no-build-isolation needs it present"
	exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

printf '== sdist offline: installing %s with cargo and pip both offline\n' "$(basename "$sdist")"
CARGO_NET_OFFLINE=true \
	python3 -m pip install --no-index --no-build-isolation \
	--target "$work/site" "$sdist"

printf '== sdist offline: the installed package computes\n'
PYTHONPATH=$work/site python3 -c "
import symfn
from symfn import s, macdonald
assert s([2, 1]) * s([1]) == s([1]) * s([2, 1])
assert macdonald.P([2]).at(q=5, t=5) == s([2]).to('m')
print('symfn', symfn.__version__, 'built and computed with the network off')
"
