#!/bin/bash
# setup_case.sh NAME MODE
#   Copy reacting case runs/rx_meml at 60 s (before its first bio update) to
#   capn-check/cases/NAME and set the biology:
#     ref  as rx_meml: MEML, production parameters, no cap_n entry
#     F    MEML with the biomass cap (set F, cap_n 5.2626)
#     G    EMM with the sharper cap (set G, cap_n 5.0932)
#   fluid_update_time 60.0001 puts the bio update on the first step after
#   60 s, as in the rx_meml smoke test. runs/rx_meml is only read.
set -uo pipefail
source /scratch/nwimer/agent/accelerate/runs/env.sh
CAP=$ROOT/capn-check
SRC=$ROOT/runs/rx_meml
NAME=${1:?name}; MODE=${2:?ref|F|G}
C=$CAP/cases/$NAME
[ -e "$C" ] && { echo "$C exists"; exit 2; }
[ -d "$SRC/processors104/60" ] || { echo "no $SRC/processors104/60"; exit 2; }
t0=$SECONDS
mkdir -p "$C/processors104" "$C/constant" || exit 2
cp -r "$SRC/0" "$SRC/system" "$SRC/Allrun" "$C/" || exit 2
# constant/ without the undecomposed mesh, which a parallel run does not
# read; Allrun only checks that constant/polyMesh is there.
for f in "$SRC"/constant/*; do
    [ "$(basename "$f")" = polyMesh ] || cp -r "$f" "$C/constant/" || exit 2
done
ln -s "$SRC/constant/polyMesh" "$C/constant/polyMesh" || exit 2
cp "$SRC"/scaler_*.txt "$SRC"/tri_*.txt "$SRC"/scripted_model_biomass_h_muconate.pt "$C/" || exit 2
cp -r "$SRC/processors104/60" "$C/processors104/" || exit 2
# The decomposed mesh is copied once and shared by the cases here.
MESH=$CAP/cases/mesh104
[ -d "$MESH" ] || cp -r "$SRC/processors104/constant" "$MESH" || exit 2
ln -s "$MESH" "$C/processors104/constant" || exit 2

D=$C/constant/microbeUpdateProperties
setv()  # setv ENTRY VALUE: replace the entry's line
{
    sed -i "s|^$1 .*|$1 $2;|" "$D"
    grep -qF "$1 $2;" "$D" || { echo "could not set $1"; exit 2; }
}
setv fluid_update_time "[ 0 0 1 0 0 0 0 ] 60.0001"
case $MODE in
ref)
    ;;
emm)  # EMM, production parameters, no cap_n entry
    setv use_ml false
    ;;
F)
    setv use_ml true
    setv y_os "[ 0 0 0 0 0 0 0 ] 2.5074"
    setv x_max "[ -1 0 0 0 1 0 0 ] 7.7983"
    setv q_max "[ -1 0 -1 0 1 0 0 ] 4.6896"
    setv K_o "[ 0 -3 0 0 1 0 0 ] 0.003125"
    setv K_s "[ 0 -3 0 0 1 0 0 ] 13.55"
    echo "cap_n           5.2626;" >> "$D"
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
    echo "unknown mode: $MODE"; exit 2
    ;;
esac
echo "$MODE biology on runs/rx_meml at 60 s, built $(date '+%F %T')" > "$C/capn_source.txt"
echo "--- $NAME ($MODE): biology entries"
grep -E "^(use_ml|y_xs|y_ms|y_os|x_max|q_max|K_o|K_s|cap_n|fluid_update_time|reaction_update_time)" "$D"
echo "--- differences from rx_meml's dictionary"
diff "$SRC/constant/microbeUpdateProperties" "$D" | grep -E '^[<>]' | cut -c1-110
echo "--- sizes"; du -sh "$C" "$MESH" | sed 's/^/    /'
echo "set up in $((SECONDS - t0)) s"
