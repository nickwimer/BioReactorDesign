#!/bin/bash
# run_case.sh NAME [END=60.002]
#   Run capn-check/cases/NAME with the cap_n build of putidaFoam on 104 tasks
#   to END (four steps; the bio update comes on the second), then print what
#   the solver read and what the bio update did.
set -uo pipefail
source /scratch/nwimer/agent/accelerate/runs/env.sh
CAP=$ROOT/capn-check
NAME=${1:?name}; END=${2:-60.002}
export PATH=$CAP/bin:$PATH
[ "$(command -v putidaFoam)" = "$CAP/bin/putidaFoam" ] \
    || { echo "wrong solver: $(command -v putidaFoam)"; exit 2; }
cd "$CAP/cases/$NAME" || exit 2
foamDictionary system/controlDict -entry endTime -set "$END" > /dev/null
t0=$SECONDS
NPROCS=104 SRUN_ARGS="-c 1" ./Allrun; rc=$?
echo "$NAME: Allrun exit $rc after $((SECONDS - t0)) s"
grep -m1 '^Exec' log.solver
grep -E "FOAM FATAL|FOAM Warning|Floating point|srun: error|Killed" log.solver | sort | uniq -c | head -5
if [ $rc -ne 0 ]; then
    grep -m1 -B2 -A14 "FOAM FATAL" log.solver
    tail -5 log.solver
fi
echo "--- bio update"
grep -m1 "Use ML for bio kinetics" log.solver
grep -m1 -A11 "MM parameters are" log.solver
grep -m1 "Avg quantities after reaction" log.solver
echo "--- timehist.dat (reaction time, O2, X, G, M before the update, then pH)"
cat timehist.dat; echo
echo "--- steps: $(grep -c '^Time = ' log.solver), last $(grep '^Time = ' log.solver | tail -1)"
echo "--- written times: $(ls processors104 | tr '\n' ' ')"
exit $rc
