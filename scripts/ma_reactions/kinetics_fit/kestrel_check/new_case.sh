#!/bin/bash
# new_case.sh NAME SRC SET
#   A standard reacting run with new kinetic parameters, for a before/after
#   against the existing runs.
#   1. runs/rx_variant.sh SRC NAME 60 3.125e-3: copy reacting case runs/SRC at
#      60 s (before its first bio update) to runs/NAME with the adopted O2
#      floor on the sink (0.1 mg/L), exactly as rx_meml_k01 was made.
#   2. Set the biology in constant/microbeUpdateProperties:
#        SET F  MEML with the biomass cap, fitted with this case's model file
#        SET G  EMM with the sharper cap
#      and fluid_update_time 60.0001, so the 4 h bio update comes on the
#      first step after 60 s and the next one would be due at 120 s.
#   The case needs the cap_n build of putidaFoam (capn-check/bin): run it
#   with capn-check/run_new_batch.sh.
set -uo pipefail
source /scratch/nwimer/agent/accelerate/runs/env.sh
NAME=${1:?name}; SRC=${2:?source case}; SET=${3:?F|G}
C=$ROOT/runs/$NAME
[ -e "$C" ] && { echo "$C exists"; exit 2; }
bash "$ROOT/runs/rx_variant.sh" "$SRC" "$NAME" 60 3.125e-3 || { echo "rx_variant.sh failed"; exit 2; }

D=$C/constant/microbeUpdateProperties
setv()  # setv ENTRY VALUE: replace the entry's line
{
    sed -i "s|^$1[ \t].*|$1 $2;|" "$D"
    grep -qF "$1 $2;" "$D" || { echo "could not set $1"; exit 2; }
}
setv fluid_update_time "[ 0 0 1 0 0 0 0 ] 60.0001"
case $SET in
F)
    setv use_ml true
    setv y_os "[ 0 0 0 0 0 0 0 ] 2.5433"
    setv x_max "[ -1 0 0 0 1 0 0 ] 7.8684"
    setv q_max "[ -1 0 -1 0 1 0 0 ] 4.7532"
    setv K_o "[ 0 -3 0 0 1 0 0 ] 0.003125"
    setv K_s "[ 0 -3 0 0 1 0 0 ] 17.527"
    echo "cap_n           3.8849;" >> "$D"
    ;;
G)
    setv use_ml false
    setv y_xs "[ 0 0 0 0 0 0 0 ] 0.031042"
    setv y_ms "[ 0 0 0 0 0 0 0 ] 0.34855"
    setv y_os "[ 0 0 0 0 0 0 0 ] 2.5168"
    setv x_max "[ -1 0 0 0 1 0 0 ] 7.7727"
    setv q_max "[ -1 0 -1 0 1 0 0 ] 4.5509"
    setv K_o "[ 0 -3 0 0 1 0 0 ] 0.003125"
    setv K_s "[ 0 -3 0 0 1 0 0 ] 13.822"
    echo "cap_n           5.0932;" >> "$D"
    ;;
*)
    echo "unknown set: $SET"; exit 2
    ;;
esac
echo "kinetics: set $SET of the kinetics fit rework (capn-check/README.txt), bio update on the first step after 60 s; needs capn-check/bin/putidaFoam; set $(date '+%F %T')" >> "$C/rx_source.txt"
echo "--- $NAME: rx_source.txt"; cat "$C/rx_source.txt"
echo "--- biology entries"
grep -E "^(use_ml|y_xs|y_ms|y_os|x_max|q_max|K_o|K_s|cap_n|fluid_update_time|reaction_update_time)" "$D"
echo "--- differences from runs/$SRC's dictionary"
diff "$ROOT/runs/$SRC/constant/microbeUpdateProperties" "$D" | grep -E '^[<>]' | cut -c1-120
echo "--- O2 sink"; grep -n -E "oxygenLimitK|useMonadSink" "$C/constant/fvModels" | head -6
echo "--- controlDict"; grep -E "^(startFrom|endTime|writeInterval|purgeWrite)" "$C/system/controlDict"
echo "--- times: $(ls "$C/processors104" | tr '\n' ' ')"
