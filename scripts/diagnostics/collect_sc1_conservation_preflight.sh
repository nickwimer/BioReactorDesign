#!/usr/bin/env bash
# Read source/case/library evidence; write only a new report in the pair root.
# This script does not build, launch CFD, alter inputs or manage allocations.

sc1_collect_evidence()
{
    git -C "$SC1_REPO" rev-parse HEAD
    cat "$SC1_AUDITROOT/source_manifest.json"
    printf '\nMATCHED-STATE MANIFEST\n'
    cat "$SC1_PAIRROOT/comparison_manifest.json"
    sha256sum "$SC1_AUDITROOT/bin/putidaFoamSC1ConservationAudit"
    for SC1_CASE in control normal; do
        printf '\nINITIAL DIRECTORIES: %s\n' "$SC1_CASE"
        ls -1 "$SC1_PAIRROOT/$SC1_CASE/processors52"
        printf '\nCHECKPOINT FIELDS: %s\n' "$SC1_CASE"
        ls -1 "$SC1_PAIRROOT/$SC1_CASE/processors52/31.199982372783058"
        printf '\nROOT ASSET HASHES: %s\n' "$SC1_CASE"
        for SC1_FILE in "$SC1_PAIRROOT/$SC1_CASE/"*; do
            if [[ -f "$SC1_FILE" ]]; then sha256sum "$SC1_FILE"; fi
        done
    done
    printf '\nCONTROL CONSTANT FILES\n'
    ls -1 "$SC1_PAIRROOT/control/constant"
    for SC1_FILE in \
        "$SC1_PAIRROOT/control/system/controlDict" \
        "$SC1_PAIRROOT/control/system/fvSchemes" \
        "$SC1_PAIRROOT/control/system/fvSolution" \
        "$SC1_PAIRROOT/control/constant/fvModels" \
        "$SC1_PAIRROOT/control/constant/fvConstraints" \
        "$SC1_PAIRROOT/control/constant/dynamicMeshDict" \
        "$SC1_PAIRROOT/control/constant/"*Properties* \
        "$SC1_PAIRROOT/control/constant/"*Transport*; do
        printf '\nFILE: %s\n' "$SC1_FILE"
        if [[ -f "$SC1_FILE" ]]; then cat "$SC1_FILE"; else echo "Not present"; fi
    done
    printf '\nCONTROL/NORMAL PHASE SETTINGS DIFFERENCES\n'
    diff -u "$SC1_PAIRROOT/control/constant/phaseProperties" "$SC1_PAIRROOT/normal/constant/phaseProperties"
    printf 'diff_status=%s (0 identical, 1 different, 2 error)\n' "$?"

    for SC1_REL in putidaFoam.C YEqns.H createFields.H bioReact.H phReact.H microbeModel.C; do
        printf '\nSOLVER SOURCE DIFF: %s\n' "$SC1_REL"
        diff -u "$SC1_REPO/OFsolvers/birdmultiphaseEulerFoam/putidaFoam/$SC1_REL" "$SC1_AUDITROOT/solver/$SC1_REL"
        printf 'diff_status=%s (0 identical, 1 different, 2 error)\n' "$?"
    done
    for SC1_REL in \
        phaseSystems/PhaseSystems/PhaseTransferPhaseSystem/PhaseTransferPhaseSystem.C \
        phaseSystems/PhaseSystems/InterfaceCompositionPhaseChangePhaseSystem/InterfaceCompositionPhaseChangePhaseSystem.C \
        phaseSystems/PhaseSystems/PopulationBalancePhaseSystem/PopulationBalancePhaseSystem.C \
        phaseSystems/phaseModel/MultiComponentPhaseModel/MultiComponentPhaseModel.C \
        phaseSystems/phaseModel/MovingPhaseModel/MovingPhaseModel.C \
        phaseSystems/phaseSystem/phaseSystemSolve.C \
        phaseSystems/diameterModels/velocityGroup/velocityGroup.C \
        phaseSystems/populationBalanceModel/populationBalanceModel/populationBalanceModel.C; do
        printf '\nDEPENDENCY SOURCE DIFF: %s\n' "$SC1_REL"
        diff -u "$SC1_REPO/OFsolvers/birdmultiphaseEulerFoam/$SC1_REL" "$SC1_AUDITROOT/$SC1_REL"
        printf 'diff_status=%s (0 identical, 1 different, 2 error)\n' "$?"
    done

    # Known trial library paths; compare them against the actual controlDict.
    for SC1_LIB in \
        /scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/libsc1HeadspaceFilmModels.so \
        /scratch/nwimer/sc1_headspace_entrainment_20261104/lib/libsc1HeadspaceEntrainmentModels.so; do
        printf '\nEXPERIMENTAL LIBRARY: %s\n' "$SC1_LIB"
        sha256sum "$SC1_LIB"
        ldd "$SC1_LIB"
    done
    echo "Evidence collection finished; review errors/differences. No scientific PASS asserted."
    return 0
}

if [[ $# -ne 2 ]]; then
    echo "Usage: bash collect_sc1_conservation_preflight.sh AUDITROOT PAIRROOT"
    false
else
    SC1_REPO=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
    SC1_AUDITROOT=$1
    SC1_PAIRROOT=$2
    SC1_REPORT="$SC1_PAIRROOT/preflight.txt"
    if [[ ! -f "$SC1_AUDITROOT/source_manifest.json" \
        || ! -x "$SC1_AUDITROOT/bin/putidaFoamSC1ConservationAudit" \
        || ! -f "$SC1_PAIRROOT/comparison_manifest.json" \
        || ! -d "$SC1_PAIRROOT/control/processors52" \
        || ! -d "$SC1_PAIRROOT/normal/processors52" \
        || -e "$SC1_REPORT" || -L "$SC1_REPORT" ]]; then
        echo "Missing build/pair evidence or report already exists; nothing overwritten."
        false
    elif (set -C; sc1_collect_evidence > "$SC1_REPORT" 2>&1); then
        printf '\nPAIRROOT=%s\nREPORT=%s\n' "$SC1_PAIRROOT" "$SC1_REPORT"
        echo "Cases prepared; evidence collected for review. No CFD launched."
        cat "$SC1_REPORT"
    else
        echo "Report creation stopped; preserve any existing report and inspect the error."
        false
    fi
fi
