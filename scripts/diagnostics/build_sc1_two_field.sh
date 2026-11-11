#!/usr/bin/env bash
# Build only. All project libraries and the nonreacting solver stay private.
sc1_two_build()
{
    if [[ "$#" -ne 1 || -z "${SLURM_JOB_ID:-}" ]]; then
        echo "Stopped: supply one fresh build path inside an active allocation."
        return 1
    fi
    local sc1_two_py="${SC1_TWO_PY:-/scratch/nwimer/.conda-envs/of-post/bin/python}"
    local sc1_two_here
    sc1_two_here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" || return 1
    if [[ ! -x "$sc1_two_py" ]]; then
        echo "Stopped: Python unavailable."
        return 1
    fi
    PYTHONDONTWRITEBYTECODE=1 "$sc1_two_py" "$sc1_two_here/sc1_two_field.py" build --out "$1"
}
sc1_two_build "$@"
