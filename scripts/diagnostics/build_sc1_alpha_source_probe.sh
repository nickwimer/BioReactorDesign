#!/usr/bin/env bash
# Run as a child with bash, not by sourcing into the allocation shell.
# Builds a fresh private observer; never prepares a case or runs CFD.

sc1_source_probe_build()
{
    if [[ "$#" -ne 1 || -z "${SLURM_JOB_ID:-}" ]]; then
        echo "Stopped: supply one fresh scratch build path inside an active allocation."
        return 1
    fi
    local sc1_probe_root="$1"
    local sc1_probe_here
    local sc1_probe_py="${SC1_PROBE_PY:-/scratch/nwimer/.conda-envs/of-post/bin/python}"
    sc1_probe_here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" || return 1
    if [[ ! -x "$sc1_probe_py" ]] || ! command -v wmake >/dev/null; then
        echo "Stopped: Python or OpenFOAM environment missing."
        return 1
    fi
    if "$sc1_probe_py" "$sc1_probe_here/sc1_alpha_source_probe.py" prepare --out "$sc1_probe_root"; then
        echo "PROBEROOT=$sc1_probe_root"
    else
        echo "Preparation stopped; preserve every existing/partial directory."
        return 1
    fi
    if (cd "$sc1_probe_root/phaseSystems" && set -C &&
        wmake libso > "$sc1_probe_root/log.build.phase" 2>&1) &&
       (cd "$sc1_probe_root/solver" && set -C &&
        wmake > "$sc1_probe_root/log.build.solver" 2>&1); then
        tail -n 12 "$sc1_probe_root/log.build.phase"
        tail -n 12 "$sc1_probe_root/log.build.solver"
        "$sc1_probe_py" "$sc1_probe_here/sc1_alpha_source_probe.py" verify --root "$sc1_probe_root"
    else
        echo "Build stopped; preserve this tree and paste the output."
        for sc1_probe_log in "$sc1_probe_root/log.build.phase" "$sc1_probe_root/log.build.solver"; do
            if [[ -f "$sc1_probe_log" ]]; then
                echo "$sc1_probe_log"
                tail -n 60 "$sc1_probe_log"
            fi
        done
        return 1
    fi
}

sc1_source_probe_build "$@"
