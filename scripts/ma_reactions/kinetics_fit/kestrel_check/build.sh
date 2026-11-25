#!/bin/bash
# build.sh: build putidaFoam with the cap_n change (capn-check/src) into
# capn-check/bin. It links the run libraries in build/two-phase-reset/lib and
# rebuilds none of them. build/two-phase-reset/bin and the workspace checkout
# are not touched.
#
# Second, compile only: the same source with wellMixed.H switched back on
# (the checkout has that include commented out), into capn-check/bin-wm, to
# check that the two changed lines of wellMixed.H compile as well.
set -uo pipefail
source /scratch/nwimer/agent/accelerate/runs/env.sh
CAP=$ROOT/capn-check
SRC=$CAP/src/birdmultiphaseEulerFoam
WM=$CAP/src-wm/birdmultiphaseEulerFoam
ORIG=$REPO/OFsolvers/birdmultiphaseEulerFoam
SIBLINGS="include phaseSystems interfacialModels interfacialCompositionModels"
N=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-8}}
mkdir -p "$CAP/bin" "$CAP/bin-wm" || exit 2
before=$(sha256sum "$BUILD/bin/putidaFoam" | cut -c1-16)

# putidaFoam's Make/options reaches its sibling directories with ../ paths:
# point them at the workspace checkout, which is only read.
for d in $SIBLINGS; do
    [ -e "$SRC/$d" ] || ln -s "$ORIG/$d" "$SRC/$d" || exit 2
done
echo "source: workspace checkout $(git -C "$REPO" log -1 --format='%h %s') + cap_n change"
echo "--- change against the workspace checkout"
diff -ru -x 'linux64*' -x lnInclude -x '*.dep' -x '*.o' "$ORIG/putidaFoam" "$SRC/putidaFoam" \
    | grep -E '^(diff|[+-][^+-])' | sed 's|^diff .*/putidaFoam/|== |' | cut -c1-140

build()  # build LABEL TREE BINDIR
{
    export FOAM_USER_APPBIN=$3
    cd "$2" || return 2
    wclean putidaFoam > /dev/null 2>&1
    local log=$CAP/log.wmake.$1 start=$SECONDS rc
    wmake -j "$N" putidaFoam > "$log" 2>&1; rc=$?
    echo "$1: wmake exit $rc after $((SECONDS-start)) s, $(grep -c -E ' error' "$log") error line(s), $(grep -ci 'warning' "$log") warning line(s)"
    echo "  compiler messages naming the change: $(grep -c -E 'cap_n|biomass_cap' "$log")"
    grep -n -E 'cap_n|biomass_cap' "$log" | head -6
    if [ $rc -ne 0 ]; then
        grep -n -m8 -E ' error' "$log"; tail -12 "$log"
    fi
    return $rc
}

echo "--- build 1: the solver as the checkout builds it, with the change"
echo "FOAM_USER_LIBBIN=$FOAM_USER_LIBBIN"
build putidaFoam "$SRC" "$CAP/bin"; rc1=$?
if [ $rc1 -eq 0 ]; then
    ls -la "$CAP/bin/putidaFoam"
    ldd "$CAP/bin/putidaFoam" | grep -E 'bird|torch|c10|not found' | sed 's/^/    /'
    "$CAP/bin/putidaFoam" -help > "$CAP/log.help" 2>&1
    echo "  -help exit $?: $(head -1 "$CAP/log.help")"
fi

echo "--- build 2 (compile only): the same with wellMixed.H switched on"
if [ ! -d "$WM/putidaFoam" ]; then
    mkdir -p "$WM" && cp -r "$SRC/putidaFoam" "$WM/putidaFoam" || exit 2
    for d in $SIBLINGS; do ln -s "$ORIG/$d" "$WM/$d" || exit 2; done
    sed -i 's|//#include "wellMixed.H"|#include "wellMixed.H"|' "$WM/putidaFoam/putidaFoam.C"
fi
grep -n 'wellMixed.H' "$WM/putidaFoam/putidaFoam.C"
build putidaFoam-wellmixed "$WM" "$CAP/bin-wm"; rc2=$?

after=$(sha256sum "$BUILD/bin/putidaFoam" | cut -c1-16)
echo "--- the solver the running cases use: $(ls -la "$BUILD/bin/putidaFoam")"
echo "    sha256 before/after this script: $before / $after"
echo "build 1 exit $rc1, build 2 exit $rc2"
exit $rc1
