#!/bin/bash
# Sets up the environment of the CI jobs (through $GITHUB_ENV and
# $GITHUB_PATH): Basilisk, the sandbox includes, Cantera (if the
# micromamba environment 'cantera' exists) and a dummy ffmpeg.

set -euo pipefail

basilisk=$HOME/basilisk/src
test -d "$basilisk" || tar xf basilisk-build.tar -C "$HOME"

# OpenSMOKE is not available: the cases use Cantera (-DUSE_OPENSMOKE=0
# in run/Makefile), but qcc still resolves the includes of the
# '#if USE_OPENSMOKE' branches, so an empty OpenSMOKE_Interface.h is
# provided. Paths are kept short: qcc crashes if the compilation command
# exceeds 1000 characters (fixed-size buffers in qcc.c). The sandbox
# headers are found through -I../src in run/Makefile.
stubs=$HOME/inc
mkdir -p "$stubs"
: > "$stubs/OpenSMOKE_Interface.h"

{
    echo "BASILISK=$basilisk"
    echo "BASILISK_INCLUDES=-I$stubs"
    echo "BASILISK_LIBS="
} >> "$GITHUB_ENV"
echo "$basilisk" >> "$GITHUB_PATH"

cantera=${MAMBA_ROOT_PREFIX:-/nonexistent}/envs/cantera
if test -f "$cantera/include/cantera/clib/ct.h"; then
    {
        echo "CANTERA_INCLUDE_PATH=$cantera/include"
        echo "CANTERA_LIBRARY_PATH=$cantera/lib"
    } >> "$GITHUB_ENV"
    grep CANTERA_VERSION "$cantera/include/cantera/base/config.h"
fi

# Movies are not checked: a dummy ffmpeg discards the frames and creates
# an empty output file (its last argument), as expected by
# $BASILISK/ppm2mp4 (a missing ffmpeg would print a warning in the log
# and break the comparison with the .ref file).
bin=$RUNNER_TEMP/bin
mkdir -p "$bin"
cat > "$bin/ffmpeg" <<'EOF'
#!/bin/sh
for a; do out=$a; done
case " $* " in *" -i - "*) cat > /dev/null ;; esac
: > "$out"
EOF
chmod +x "$bin/ffmpeg"
echo "$bin" >> "$GITHUB_PATH"
