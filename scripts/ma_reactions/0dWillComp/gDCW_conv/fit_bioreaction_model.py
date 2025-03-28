import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize


def setmodelconstants(p):

    modelConstants = {}
    # modelConstants["Fs_max"] = p[0]  # molS/m^3 / kgBio/m^3
    modelConstants["q_max"] = p[0]  # molS/m^3 / kgBio/m^3
    modelConstants["bio_max"] = 7.9  # kg/m^3

    modelConstants["K_o"] = p[1]  # mol/m^3
    # modelConstants['K_s'] = p[3] # mol/m^3

    modelConstants["Y_xs"] = p[2]  # g/molS
    modelConstants["Y_ms"] = p[3]  # molB/molS
    modelConstants["Y_os"] = p[4]  # molO/molS

    modelConstants["o2sat"] = 0.214
    modelConstants["K_s"] = p[5] # 0.92  # mol/m^3

    modelConstants["k_la"] = 50 # p[3]
    return modelConstants


def ode_model(y, t, mc):
    """
    Definition of simplified metabolic model, as defined
    in the document.

    Inputs:
    --f is the concentrations of the various species:
        bio, glu, muc = f
        all in mol/m^3, except biomass, which is kg/m^3
    --otr is the oxygen transfer rate, in mM/h
    --mc is the set of model constants, in a dict
    Outputs:
    --dfdt is the model output
    """
    bio, glu, muc, o2 = y

    # fixed k_la
    k_la = mc["k_la"]

    # Extract model constants
    F_s = glu / (glu + mc["K_s"])
    F_o = o2 / (o2 + mc["K_o"])

    # Calculate qs
    qs = mc["q_max"]*F_s * F_o

    # Calculate rates
    rbio = mc["Y_xs"] * qs * bio * (1 - bio / mc["bio_max"])  # biomass
    rglu = -qs * bio  # glucose
    rmuc = mc["Y_ms"] * qs * bio  # muconate

    otr = k_la * (mc["o2sat"] - o2)  # oxygen transfer rate
    roxy = -mc["Y_os"] * qs * bio + otr  # oxygen

    return np.array([rbio, rglu, rmuc, roxy])


def loss_function(params, init_cond, t_final, exp_data, normalize=False):
    """Compute the mse between the model and experimental data"""

    mc = setmodelconstants(params)

    # Just solve at the experimental time points
    time_int = exp_data["time"]

    sol_y = odeint(ode_model, init_cond, time_int, args=(mc,))

    F_s = sol_y[:, 1] / (sol_y[:, 1] + mc["K_s"])
    F_o = sol_y[:, 3] / (sol_y[:, 3] + mc["K_o"])

    qs = mc["q_max"] * F_s * F_o

    our = qs * mc["Y_os"] * sol_y[:, 0]

    # Compute the mean squared error
    if normalize:
        mse = np.mean((exp_data["bio"] - sol_y[:, 0]) ** 2) / np.mean(
            exp_data["bio"] ** 2
        )
        mse += np.mean((exp_data["glu"] - sol_y[:, 1]) ** 2) / np.mean(
            exp_data["glu"] ** 2
        )
        mse += np.mean((exp_data["muc"] - sol_y[:, 2]) ** 2) / np.mean(
            exp_data["muc"] ** 2
        )
        mse += np.mean((exp_data["o2"] - sol_y[:, 3]) ** 2) / np.mean(
            exp_data["o2"] ** 2
         )
    else:
        mse = np.mean((exp_data["bio"] - sol_y[:, 0]) ** 2)
        mse += np.mean((exp_data["glu"] - sol_y[:, 1]) ** 2)
        mse += np.mean((exp_data["muc"] - sol_y[:, 2]) ** 2)
        mse += np.mean((exp_data["o2"] - sol_y[:, 3]) ** 2)

    return mse

def integrate_solution(params, init_cond, t_final, exp_data):
    """Integrate the ODE model with the given parameters"""

    # Make a smooth time integration
    time_int = np.linspace(0, t_final, 1000)

    mc = setmodelconstants(params)
    print(mc)

    sol_y = odeint(ode_model, init_cond, time_int, args=(mc,))
    print(sol_y)

    F_s = sol_y[:, 1] / (sol_y[:, 1] + mc["K_s"])
    F_o = sol_y[:, 3] / (sol_y[:, 3] + mc["K_o"])

    qs = mc["q_max"] * F_s * F_o

    our = qs * mc["Y_os"] * sol_y[:, 0]

    # Save results that correspond to experimental data points
    sol_y_int = np.array([np.interp(exp_data["time"], time_int, y) for y in np.transpose(sol_y)])
    np.savetxt(os.path.join(case_dir,"soln.dat"), np.transpose(np.vstack((exp_data["time"], sol_y_int[:-1]))), delimiter=" ")

    # Plot results
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), squeeze=False)
    axs[0, 0].plot(exp_data["time"], exp_data["bio"], "ro", label="expr")
    axs[0, 0].plot(time_int, sol_y[:, 0], "b-", label="model")
    axs[0, 0].set_ylabel(r"biomass ($kg/m^3$)", fontsize=14)
    axs[0, 0].set_xlabel("time (h)", fontsize=14)

    axs[0, 1].plot(exp_data["time"], exp_data["glu"], "ro", label="expr")
    axs[0, 1].plot(time_int, sol_y[:, 1], "b-", label="model")
    axs[0, 1].set_ylabel(r"glucose ($mol/m^3$)", fontsize=14)
    axs[0, 1].set_xlabel("time (h)", fontsize=14)

    axs[1, 0].plot(exp_data["time"], exp_data["muc"], "ro", label="expr")
    axs[1, 0].plot(time_int, sol_y[:, 2], "b-", label="model")
    axs[1, 0].set_ylabel(r"muconate ($mol/m^3$)", fontsize=14)
    axs[1, 0].set_xlabel("time (h)", fontsize=14)

    axs[1, 1].plot(exp_data["time"], exp_data["o2"], "ro", label="expr")
    axs[1, 1].plot(time_int, sol_y[:, 3], "b-", label="model")
    axs[1, 1].set_ylabel("dissolved O2 ($mol/m^3$)", fontsize=14)
    axs[1, 1].set_xlabel("time (h)", fontsize=14)

    return fig

def plot_rates(params, init_cond, t_final):
    """Integrate the ODE model with the given parameters"""

    rates = np.empty((1000,4))
    
    # Make a smooth time integration
    time_int = np.linspace(0, t_final, 1000)

    mc = setmodelconstants(params)
    print(mc)

    sol_y = odeint(ode_model, init_cond, time_int, args=(mc,))
    print(sol_y)
    print(sol_y.shape)
    print(rates.shape)
    

    for i in range(1000):
        rates[i] = ode_model(sol_y[i], time_int, mc)

    F_s = sol_y[:, 1] / (sol_y[:, 1] + mc["K_s"])
    F_o = sol_y[:, 3] / (sol_y[:, 3] + mc["K_o"])

    qs = mc["q_max"] * F_s * F_o

    our = qs * mc["Y_os"] * sol_y[:, 0]
    
    # Plot results
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), squeeze=False)
    axs[0, 0].plot(time_int, rates[:, 0]/sol_y[:, 0], "b-", label="model")
    axs[0, 0].set_ylabel(r"biomass ($1/h$)", fontsize=14)
    axs[0, 0].set_xlabel("time (h)", fontsize=14)

    axs[0, 1].plot(time_int, rates[:, 1]/sol_y[:, 0], "b-", label="model")
    axs[0, 1].set_ylabel(r"glucose ($mmol/gDCW/h$)", fontsize=14)
    axs[0, 1].set_xlabel("time (h)", fontsize=14)

    axs[1, 0].plot(time_int, rates[:, 2]/sol_y[:, 0], "b-", label="model")
    axs[1, 0].set_ylabel(r"muconate ($mmol/gDCW/h$)", fontsize=14)
    axs[1, 0].set_xlabel("time (h)", fontsize=14)

    axs[1, 1].plot(time_int, -our/sol_y[:, 0], "b-", label="model")
    axs[1, 1].set_ylabel("-OUR ($mmol/gDCW/h$)", fontsize=14)
    axs[1, 1].set_xlabel("time (h)", fontsize=14)

    return fig

if __name__ == "__main__":

    # Case directory
    parser = argparse.ArgumentParser(description="0d bioreaction kinetic model of P. putida")
    parser.add_argument(
        "-d",
        "--directory",
        dest="fdir",
        help="directory containing experiment data, initial conditions, and initial parameter guess",
        type=str,
        default=".",
    )
    parser.add_argument(
        "-pf",
        "--param_fit",
        action=argparse.BooleanOptionalAction,
        help="Perform parameter fitting based on initial guess provided?",
    )
    args = parser.parse_args()
    case_dir = os.path.abspath(args.fdir)

    # Data order: t, bio, glu, muc, o2
    exp_data = np.loadtxt(os.path.join(case_dir, "exptdata_ma.dat"))
    init_cond = np.loadtxt(os.path.join(case_dir, "initial_cond.dat"))

    # Extract experimental data
    exp_data_dict = {
        "time": exp_data[:, 0],
        "bio": exp_data[:, 1],
        "glu": exp_data[:, 2],
        "muc": exp_data[:, 3],
        "o2": exp_data[:, 4],
    }
    # print(exp_data_dict)

    # Initial parameter order: q_max, K_o, Y_xs, Y_ms, Y_os
    p0 = np.loadtxt(os.path.join(case_dir, "initial_param.dat"))

    if args.param_fit:
        # Solve ODE
        fig = integrate_solution(p0, init_cond, 33, exp_data_dict)
        fig.suptitle(f"Initial Parameters: {p0}", fontsize=16)
        fig.savefig(os.path.join(case_dir, "initial_parameters.png"))
        plt.close(fig)

        # Optimize parameters
        initial_guess = p0
        result = minimize(
            loss_function,
            initial_guess,
            bounds=[(0, None)] * 6,
            args=(init_cond, 33, exp_data_dict, True),  # Final argument is for norm mse
        )
        fitted_params = result.x

        # Print results
        print("Initial Parameters:", p0)
        print("Fitted Parameters:", fitted_params)
    else:
        fitted_params = p0

    np.savetxt(os.path.join(case_dir,"fitted_param.dat"), np.transpose(fitted_params), delimiter=" ")

    # Solve ODE with fitted parameters
    fig = integrate_solution(fitted_params, init_cond, 33, exp_data_dict)
    fig.suptitle(f"Fitted Parameters: {fitted_params}", fontsize=16)
    fig.savefig(os.path.join(case_dir, "fitted_parameters.png"))
    plt.close(fig)

    plt.show()

    fig = plot_rates(fitted_params, init_cond, 33)
    fig.suptitle("Rates", fontsize=16)
    fig.savefig(os.path.join(case_dir, "CFD_rates.png"))
    plt.close(fig)

    
