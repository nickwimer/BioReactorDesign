#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=104
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/nwimer/agent/accelerate/runs/%x-%j.out
# run_new_batch.sh NAME END   run (or continue) runs/NAME to END on one node,
# as runs/run_case_batch.sh does, but with the cap_n build of putidaFoam
# (capn-check/bin). The case's Allrun takes the solver from BIRD_SOLVER.
set -uo pipefail
export BIRD_SOLVER=/scratch/nwimer/agent/accelerate/capn-check/bin/putidaFoam
[ -x "$BIRD_SOLVER" ] || { echo "no solver at $BIRD_SOLVER"; exit 2; }
"/scratch/nwimer/agent/accelerate/runs/stage1.sh" "${1:?name}" run "${2:?end time}" 104
