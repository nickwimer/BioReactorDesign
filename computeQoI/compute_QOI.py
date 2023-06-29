import argparse
import sys

import numpy as np

sys.path.append("util")
import os
import pickle

from bcr_util import *
from ofio import *

parser = argparse.ArgumentParser(
    description="Compute means QoI of OpenFOAM fields"
)
parser.add_argument(
    "-f",
    "--caseFolder",
    type=str,
    metavar="",
    required=True,
    help="caseFolder to analyze",
    default=None,
)
parser.add_argument(
    "-vl",
    "--var_list",
    nargs="+",
    help="List of variables to compute",
    default=[
        "GH",
        "GH_height",
        "d",
        "CO2_liq",
        "CO_liq",
        "H2_liq",
        "kla_CO2",
        "kla_CO",
        "kla_H2",
    ],
    required=False,
)
parser.add_argument(
    "-avg",
    "--windowAve",
    type=int,
    metavar="",
    required=False,
    help="Window Average",
    default=1,
)
parser.add_argument(
    "-conv",
    "--windowConv",
    type=int,
    metavar="",
    required=False,
    help="Window Convergence",
    default=1,
)
args = parser.parse_args()


case_path = args.caseFolder
var_name_list = args.var_list
time_float_sorted, time_str_sorted = getCaseTimes(case_path)
mesh_time_str = getMeshTime(case_path)
cellCentres = readMesh(
    os.path.join(case_path, f"meshCellCentres_{mesh_time_str}.obj")
)
nCells = len(cellCentres)

window_ave = min(args.windowAve, len(time_str_sorted) - 1)
window_conv = min(args.windowConv, len(time_str_sorted) - 1)

variables = {}
variables_conv = {}
for name in var_name_list:
    variables_conv[name] = {}
    variables_conv[name]["x"] = []
    variables_conv[name]["y"] = []


def get_var(
    case_path, time_folder, mesh_time_str, cellCentres, nCells, val_dict, name
):
    localFolder = os.path.join(case_path, time_folder)
    localFolder_vol = os.path.join(case_path, mesh_time_str)
    if name == "GH":
        var, val_dict = computeGH(
            localFolder, localFolder_vol, nCells, cellCentres, val_dict
        )
    elif name == "GH_height":
        var, val_dict = computeGH_height(
            localFolder,
            nCells,
            cellCentres,
            height_liq_base=7.0,
            val_dict=val_dict,
        )
    elif name == "d":
        var, val_dict = computeDiam(localFolder, nCells, cellCentres, val_dict)
    elif name == "CO2_liq":
        var, val_dict = computeSpec_liq(
            localFolder,
            nCells,
            field_name="CO2.liquid",
            key="co2_liq",
            cellCentres=cellCentres,
            val_dict=val_dict,
        )
    elif name == "CO_liq":
        var, val_dict = computeSpec_liq(
            localFolder,
            nCells,
            field_name="CO.liquid",
            key="co_liq",
            cellCentres=cellCentres,
            val_dict=val_dict,
        )
    elif name == "H2_liq":
        var, val_dict = computeSpec_liq(
            localFolder,
            nCells,
            field_name="H2.liquid",
            key="h2_liq",
            cellCentres=cellCentres,
            val_dict=val_dict,
        )
    elif name == "kla_CO":
        var, val_dict = computeSpec_kla(
            localFolder,
            nCells,
            key_suffix="co",
            cellCentres=cellCentres,
            val_dict=val_dict,
        )
    elif name == "kla_CO2":
        var, val_dict = computeSpec_kla(
            localFolder,
            nCells,
            key_suffix="co2",
            cellCentres=cellCentres,
            val_dict=val_dict,
        )
    elif name == "kla_H2":
        var, val_dict = computeSpec_kla(
            localFolder,
            nCells,
            key_suffix="h2",
            cellCentres=cellCentres,
            val_dict=val_dict,
        )
    else:
        sys.exit(f"ERROR: unknown variable {name}")

    return var, val_dict


print(f"Case : {case_path}")

print(f"Window Ave")
for i_ave in range(window_ave):
    time_folder = time_str_sorted[-i_ave - 1]
    print(f"\tTime : {time_folder}")
    case_variables = []
    val_dict = {}
    for name in var_name_list:
        var, val_dict = get_var(
            case_path,
            time_folder,
            mesh_time_str,
            cellCentres,
            nCells,
            val_dict=val_dict,
            name=name,
        )

        if i_ave == 0:
            variables[name] = var / window_ave
        else:
            variables[name] += var / window_ave


print(f"Window Conv")
for i_conv in range(window_conv):
    time_folder = time_str_sorted[-window_conv + i_conv]
    print(f"\tTime : {time_folder}")
    case_variables = []
    val_dict = {}
    for name in var_name_list:
        var, val_dict = get_var(
            case_path,
            time_folder,
            mesh_time_str,
            cellCentres,
            nCells,
            val_dict=val_dict,
            name=name,
        )
        variables_conv[name]["x"] += [float(time_folder)]
        variables_conv[name]["y"] += [var]


with open(os.path.join(case_path, "qoi.pkl"), "wb") as f:
    pickle.dump(variables, f)
with open(os.path.join(case_path, "qoi_conv.pkl"), "wb") as f:
    pickle.dump(variables_conv, f)
