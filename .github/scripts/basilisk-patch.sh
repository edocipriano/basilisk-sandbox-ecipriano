#!/bin/bash
# Generates .github/basilisk.patch, applied by setup-env.sh to the official
# Basilisk in the CI jobs, so that the tests use the same headers as the
# local Basilisk.
#
# Usage: .github/scripts/basilisk-patch.sh [local Basilisk src, default $BASILISK]
#
# The patch is a unified diff (patch -p1 -d $BASILISK) of the headers which
# exist in the official version and differ in the local one (test/,
# examples/ and darcsit/ are ignored). Headers which exist only locally are
# not included. Run it again after recording a change in the local Basilisk,
# or when the CI fails because upstream already contains a change.

set -euo pipefail

local=$(cd "${1:-$BASILISK}" && pwd)
url=${BASILISK_URL:-https://basilisk.fr/basilisk/basilisk.tar.gz}
out=$(cd "$(dirname "$0")/.." && pwd)/basilisk.patch

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
curl -sSfL --connect-timeout 15 --retry 3 "$url" | tar xz -C "$tmp"
official=$tmp/basilisk/src

cd "$official"
: > "$tmp/patch"
find . -name '*.h' \
     -not -path './test/*' -not -path './examples/*' -not -path './darcsit/*' \
    | sed 's|^\./||' | LC_ALL=C sort | while read -r f; do
    test -f "$local/$f" || continue
    if ! cmp -s "$f" "$local/$f"; then
        diff -u --label "a/$f" --label "b/$f" "$f" "$local/$f" \
             >> "$tmp/patch" || true
        echo "$f" >&2
    fi
done

if test -s "$tmp/patch"; then
    mv "$tmp/patch" "$out"
    echo "written $out" >&2
else
    rm -f "$out"
    echo "local headers identical to the official ones: no patch" >&2
fi
