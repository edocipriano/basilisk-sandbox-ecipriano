#!/bin/bash
# Selects the cases in run/ (with a .ref file) which need to be run.
#
# As in Basilisk's runtest (which skips a test if the checksum of the
# compiled <case>.s did not change), a case is skipped if it already
# passed with exactly the same inputs: the key of a case is a hash of
# the contents of <case>.c, <case>.ref, all the headers it includes (as
# listed by 'make <case>.s.d', i.e. 'qcc -MD'), run/Makefile, the CI
# scripts, the Cantera packages, the runner image and the Basilisk
# version. Passed keys are stored as (tiny) GitHub cache entries, so
# that this works across branches.
#
# Environment: BASILISK, BASILISK_INCLUDES, CANTERA_*, STAMP (Basilisk
# version), CANTERA (conda packages), RUNNER (runner image), FORCE
# (non-empty to run everything), GH_TOKEN, GH_REPO.
#
# Prints, for $GITHUB_OUTPUT, the JSON matrices of the cases to run and
# of those which are up to date:
#   cases=[{"case": ..., "key": ..., "cantera": ...}, ...]
#   uptodate=[{"case": ..., "key": ...}, ...]

set -euo pipefail
cd "$(dirname "$0")/../../run"

common=$( (cat Makefile ../.github/scripts/select-cases.sh ../.github/scripts/setup-env.sh
           echo "$STAMP ${CANTERA:-} ${RUNNER:-}") | sha256sum)

passed=$(gh cache list --key pass- --limit 1000 --json key --jq '.[].key' || true)

torun=()
uptodate=()
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
        uptodate+=("{\"case\":\"$c\",\"key\":\"$key\"}")
    else
        echo "to run: $c ($key, $(wc -w <<< "$deps") dependencies, cantera=$cantera)" >&2
        torun+=("{\"case\":\"$c\",\"key\":\"$key\",\"cantera\":$cantera}")
    fi
done

(IFS=,; echo "cases=[${torun[*]:-}]"; echo "uptodate=[${uptodate[*]:-}]")
