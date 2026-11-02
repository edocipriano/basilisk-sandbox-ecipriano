#!/bin/bash
# Selects the cases in run/ (with a .ref file) which need to be run.
#
# As in Basilisk's runtest (which skips a test if the checksum of the
# compiled <case>.s did not change), a case is skipped if it already
# passed with exactly the same inputs: the key of a case is a hash of
# the contents of <case>.c, <case>.ref, all the headers it includes (as
# listed by 'make <case>.s.d', i.e. 'qcc -MD'), run/Makefile, the
# workflow and its scripts (which fix the Cantera version) and the
# Basilisk version. Passed keys are stored as (tiny) GitHub cache
# entries, so that this works across branches.
#
# Environment: BASILISK, BASILISK_INCLUDES, CANTERA_*, STAMP (Basilisk
# version), FORCE (non-empty to run everything), GH_TOKEN, GH_REPO.
#
# Prints the JSON matrix [{"case": ..., "key": ..., "cantera": ...}, ...]
# on stdout.

set -euo pipefail
cd "$(dirname "$0")/../../run"

common=$( (cat Makefile ../.github/workflows/tests.yml \
               ../.github/scripts/select-cases.sh ../.github/scripts/setup-env.sh
            echo "$STAMP") | sha256sum)

passed=$(gh cache list --key pass- --limit 1000 --json key --jq '.[].key' || true)

items=()
for ref in *.ref; do
    c=${ref%.ref}
    cantera=false
    if grep -qE '^\s*#\s*include\s*"cantera/' "$c.c"; then
        cantera=true
    fi
    rm -f "$c.s.d"
    make -s "$c.s.d" >&2 || echo "warning: could not list dependencies of $c" >&2
    deps=$(tr ' \\' '\n\n' < "$c.s.d" 2> /dev/null | grep -v -e '^$' -e ':$' | sort -u || true)
    hash=$( (echo "$common"
             sha256sum < "$c.c"
             sha256sum < "$c.ref"
             for f in $deps; do
                 if test -f "$f"; then sha256sum < "$f"; fi
             done) | sha256sum | cut -c1-16)
    key="pass-$c-$hash"
    if test -z "${FORCE:-}" && grep -qx "$key" <<< "$passed"; then
        echo "up to date: $c ($key)" >&2
    else
        echo "to run: $c ($key, $(wc -w <<< "$deps") dependencies, cantera=$cantera)" >&2
        items+=("{\"case\":\"$c\",\"key\":\"$key\",\"cantera\":$cantera}")
    fi
done

(IFS=,; echo "[${items[*]:-}]")
