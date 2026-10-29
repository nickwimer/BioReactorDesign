#!/usr/bin/env bash
# Prepare a diagnostic-only SC1 continuation without touching production files.
# Usage: bash prepare_sc1_bubble_audit.sh LOCAL_SOLVER_SRC PRODUCTION_CASE MODEL_BASE SCRATCH_DEST
# Requires an initialized OpenFOAM 9 environment with foamDictionary and wmake.
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: bash $0 LOCAL_SOLVER_SRC PRODUCTION_CASE MODEL_BASE SCRATCH_DEST" >&2
    exit 2
fi

SRC=$(readlink -f "$1")
CASE=$(readlink -f "$2")
BASE=$(readlink -f "$3")
DEST="$4"

if [[ "$DEST" != /scratch/nwimer/* || -e "$DEST" || -L "$DEST" ]]; then
    echo "Refusing destination: must be a NEW directory below /scratch/nwimer/" >&2
    exit 2
fi

command -v foamDictionary >/dev/null || {
    echo "First source your RHEL9 OpenFOAM environment (foamDictionary missing)." >&2
    exit 2
}

for f in \
    "$SRC/putidaFoam.C" "$SRC/Make/files" \
    "$CASE/system/controlDict" "$CASE/system/fvSolution" \
    "$CASE/constant/phaseProperties" \
    "$CASE/processors52/constant" "$CASE/processors52/30"
do
    if [[ ! -e "$f" ]]; then
        echo "Missing required source/case path: $f" >&2
        exit 2
    fi
done

for a in scripted_model_biomass_h_muconate.pt \
         scaler_x_mean.txt scaler_x_scale.txt scaler_y_mean.txt scaler_y_scale.txt \
         tri_points.txt tri_simplices.txt
do
    if [[ ! -f "$CASE/$a" && ! -f "$BASE/$a" ]]; then
        echo "Missing required model asset in both CASE and BASE: $a" >&2
        exit 2
    fi
done

REPO=$(cd "$(dirname "$0")/../.." && pwd)
HEADER="$REPO/OFsolvers/birdmultiphaseEulerFoam/putidaFoam/bubbleMassAudit.H"
test -f "$HEADER" || { echo "Missing diagnostic header $HEADER" >&2; exit 2; }

mkdir -p "$DEST"
mkdir -p "$DEST/run/processors52"
cp -a "$SRC" "$DEST/solver"
cp -a "$CASE/system" "$CASE/constant" "$DEST/run/"
cp -a "$CASE/processors52/constant" "$CASE/processors52/30" "$DEST/run/processors52/"
cp "$HEADER" "$DEST/solver/bubbleMassAudit.H"

for a in scripted_model_biomass_h_muconate.pt \
         scaler_x_mean.txt scaler_x_scale.txt scaler_y_mean.txt scaler_y_scale.txt \
         tri_points.txt tri_simplices.txt
do
    origin="$BASE/$a"
    [[ ! -f "$CASE/$a" ]] || origin="$CASE/$a"
    cp -p "$origin" "$DEST/run/$a"
done

python3 - "$DEST/solver/putidaFoam.C" "$DEST/solver/Make/files" <<'PY'
from pathlib import Path
import re
import sys

solver = Path(sys.argv[1])
make = Path(sys.argv[2])
src = solver.read_text()

needle = "        runTime.write();"
if src.count(needle) != 1 or "bubbleMassAudit.H" in src:
    raise SystemExit(
        "Cannot safely patch source: runTime.write() must occur once, "
        "and the diagnostic must not already be present."
    )
hook = """        // Opt-in, write-time-only bubble and headspace mass audit.
        if
        (
            pimple.dict().lookupOrDefault<Switch>
            (
                "bubbleMassAudit",
                false
            )
         && runTime.writeTime()
        )
        {
            #include "bubbleMassAudit.H"
        }

"""
solver.write_text(src.replace(needle, hook + needle))

cfg = make.read_text()
new_cfg, n = re.subn(
    r"(?m)^[ \t]*EXE[ \t]*=.*$",
    "EXE = $(FOAM_USER_APPBIN)/putidaFoamMask01BubbleAudit",
    cfg,
)
if n != 1:
    raise SystemExit("Cannot safely set unique Make/files EXE output")
make.write_text(new_cfg)
PY

RUN="$DEST/run"
foamDictionary "$RUN/system/controlDict" -entry startFrom -set latestTime
foamDictionary "$RUN/system/controlDict" -entry startTime -set 30
foamDictionary "$RUN/system/controlDict" -entry endTime -set 30.3
foamDictionary "$RUN/system/controlDict" -entry writeControl -set runTime
foamDictionary "$RUN/system/controlDict" -entry writeInterval -set 0.1
foamDictionary "$RUN/system/fvSolution" -entry PIMPLE/bubbleMassAudit -set true

echo
echo "Diagnostic copy prepared:"
echo "  Source: $DEST/solver"
echo "  Case:   $DEST/run"
echo "  Build:  (cd '$DEST/solver' && wmake)"
echo "  Run:    (cd '$DEST/run' && srun -n 52 putidaFoamMask01BubbleAudit -parallel -fileHandler collated > log.bubbleAudit 2>&1)"
echo "  Inspect: grep -E 'BUBBLE_MASS_AUDIT|HEADSPACE_PRESSURE_AUDIT|FOAM FATAL|Error' '$DEST/run/log.bubbleAudit'"
echo
echo "BEFORE RUNNING: inspect cloned controlDict, model/scaler hashes, "
echo "and local RHEL9 solver divergence from the GitHub diagnostic branch."
