import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate as igt
from scipy import optimize as opt
from scipy.special import gammainc


O2IND = 0    
XIND = 1
GLUCIND = 2
MAIND = 3
NVAR = MAIND + 1
varnames = ["oxygen", "biomass", "glucose", "muconate"]


def dfdt(f, t, kLa, mc):
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
    o2, bio, glu, muc = f

    F_s = mc["Fs_max"] * glu / (glu + mc["K_s"])
    F_o = mc["Fo_max"] * o2  / (o2 + mc["K_o"])

    qs = F_s * F_o

    rbio = mc["Y_xs"] * qs * bio * (1 - bio / mc["bio_max"])
    rg = -qs * bio
    rm = mc["Y_ms"] * qs * bio

    otr = kLa * (mc["o2sat"] - o2)
    ro = -mc["Y_os"] * qs * bio + otr

    return np.array([ro, rbio, rg, rm])


def setmodelconstants(p):

    modelConstants = {}
    modelConstants["Fs_max"] = p[0] # molS/m^3 / kgBio/m^3
    modelConstants["Fo_max"] = p[1] # molS/m^3 / kgBio/m^3
    modelConstants["bio_max"] = p[2] # kg/m^3

    modelConstants["K_o"] = p[3]  # mol/m^3
    # modelConstants['K_s'] = p[3] # mol/m^3

    modelConstants["Y_xs"] = p[4]  # g/molS
    modelConstants["Y_ms"] = p[5]  # molB/molS
    modelConstants["Y_os"] = p[6]  # molO/molS

    modelConstants["o2sat"] = 0.214
    modelConstants["K_s"] = 0.92  # mol/m^3

    return modelConstants


def getrmserror(f, t, exptdata, ind1, ind2):

    # bdo data
    t_dat = exptdata[:, ind2 - 1]
    exp_dat = exptdata[:, ind2]
    rmserr = 0.0
    interp_data = np.interp(t_dat, t, f[:, ind1])
    rmserr = np.sqrt(np.mean((interp_data - exp_dat) ** 2))
    return rmserr


def residual(p, f0, kLa, tfinal, exptdata):

    t = np.linspace(0, tfinal, 200)
    mc = setmodelconstants(p)
    f = igt.odeint(dfdt, f0, t, args=(kLa, mc))

    globalrmserr = 0.0
    weights = np.zeros(NVAR) + 1.0
    for i in range(1, NVAR):
        globalrmserr += weights[i] * getrmserror(f, t, exptdata, i, 2 * i - 1)
    # print(globalrmserr)
    # rmserr=getrmserror(f,t,exptdata,BDOIND,BDOIND_DAT)
    # print(rmserr)
    return globalrmserr

    # np.savetxt("soln.dat",np.transpose(np.vstack((t,np.transpose(f)))),delimiter=" ")


def integratesoln(p, f0, kLa, tfinal, exptdata, avg_comp, cell_comp, ax, nrows, ncols):

    t = np.linspace(0, tfinal, 1000)
    mc = setmodelconstants(p)
    f = igt.odeint(dfdt, f0, t, args=(kLa, mc))

    globalrmserr = 0.0
    weights = np.zeros(NVAR) + 1.0
    for i in range(1, NVAR):
        err = weights[i] * getrmserror(f, t, exptdata, i, 2 * i - 1)
        print("solnerr for %d=%e" % (i, err))
        globalrmserr += err

    print("globalerr:", globalrmserr)
    np.savetxt("soln.dat", np.transpose(np.vstack((t, np.transpose(f)))), delimiter=" ")
    for i in range(1, NVAR):
        # plt.figure()
        row = int((i - 1) / (ncols-2))
        col = int((i - 1) % ncols)
        
        ax[0][i-1].set_ylabel("["+varnames[i]+"] mol/m^3")
        ax[0][i-1].set_xlabel("time (h)")
        ax[0][i-1].plot(exptdata[:, 2 * i - 2], exptdata[:, 2 * i - 1], "ro", label="expr")
        ax[0][i-1].plot(t, f[:, i],"gv", fillstyle="none", label="0d")
        ax[0][i-1].plot(avg_comp[:, 0], avg_comp[:, i+1], "k*", label="avged chem")
        ax[0][i-1].plot(cell_comp[:, 0], cell_comp[:, i+1], color="tab:orange", label="cell-wise chem")
        if(i == 1):
            fig.legend(loc="upper right")


if __name__ == "__main__":

    # initial condition order: O2, bio, glu, muc
    f0 = np.array([0.214, 0.5, 73.0, 0])
    kLa = 50
    tfinal = 30

    # parameter order: t, bio, t, glu, t, muc, t, our
    exptdata = np.loadtxt("exptdata_ma.csv")
    avg_comp = np.loadtxt("timehist_avg.dat")
    cell_comp = np.loadtxt("timehist_cell.dat")
    # parameter order: Fs_max, Fo_max, bio_max, K_o, Y_xs, Y_ms, Y_os
    # knobs that can be fiddled with: Y_os, Fo_max, Fs_max, K_o, kLa
    p0 = np.array([1.1, 1, 7.9, 0.0214, 0.109, 0.3, 0.0467])

    nrows = 1
    ncols = 3
    fig = plt.figure(figsize=(15, 5))
    ax = fig.subplots(1, 3, squeeze=False)
    # (fig, ax) = plt.subplots(nrows, ncols)
    integratesoln(p0, f0, kLa, tfinal, exptdata, avg_comp, cell_comp,  ax, nrows, ncols)
    fig.suptitle('Well-mixed P. putida reaction advance for different chemistry implementations')
    plt.savefig('chem_comp.pdf')
    plt.show()
