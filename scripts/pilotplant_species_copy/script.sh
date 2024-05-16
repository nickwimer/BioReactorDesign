#!/bin/bash
#SBATCH --job-name=Biomethan
##SBATCH --partition=standard
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=50  # use 52
#SBATCH --time=2-00:00:00
#SBATCH --account=hdcomb
#SBATCH --output=log.out
#SBATCH --error=log.err

module purge
ml PrgEnv-cray
source /projects/gas2fuels/ofoam_cray_mpich/OpenFOAM-dev/etc/bashrc

. ./presteps.sh
decomposePar -latestTime -fileHandler collated 
srun -n 400 multiphaseEulerFoam -parallel -fileHandler collated

# post-process
reconstructPar -newTimes
