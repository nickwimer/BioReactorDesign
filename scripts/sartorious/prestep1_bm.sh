#!/bin/bash
# source activate base
# conda activate bird
rm ./blockMeshDict_reactor
rm -r dynamicCode
rm -r 0
cp -r 0.org 0
# python3 system/write_bmesh_file.py
# conda deactivate
# source /scratch/jream/OpenFOAM/OpenFOAM-9/etc/bashrc
# of9
# cd run/putida_stir_small/pH_test
# blockMesh -dict ./blockMeshDict_reactor
# stitchMesh -perfect -overwrite inside_to_hub inside_to_hub_copy
# stitchMesh -perfect -overwrite hub_to_rotor hub_to_rotor_copy
# transformPoints "rotate=((0 0 1)(0 1 0))"
# surfaceToPatch -tol 1e-3 sparger.stl
# export newmeshdir=$(foamListTimes -latestTime)
# rm -rf constant/polyMesh/
# cp -r $newmeshdir/polyMesh ./constant
# rm -rf $newmeshdir
# sed -i '' -e 's/patch0/inlet/g' ./constant/polyMesh/boundary
# sed -i '' -e 's/zone0/inlet/g' ./constant/polyMesh/boundary
# rm -rf 0
# cp -r 0.org 0
# setFields
