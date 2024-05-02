import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate as igt
from scipy import optimize as opt
from scipy.special import gammainc


O2IND = 0    # leaving o2 in for now as electron sink (F_e); don't know if accurate
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

    F_s = glu / (glu + mc["K_s"])
    F_e = o2  / (o2 + mc["K_e"])

    qs = mc["qsmx"] * F_s * F_e

    rbio = mc["Y_xs"] * qs * bio * (1 - bio / mc["bio_mx"])
    rg = -qs * bio
    rm = mc["Y_ms"] * qs * bio

    otr = kLa * (mc["o2sat"] - o2)
    ro = -mc["Y_os"] * qs * bio + otr

    return np.array([ro, rbio, rg, rm])


def setmodelconstants(p):

    modelConstants = {}
    modelConstants["qsmx"] = p[0]  # molS/m^3 / kgBio/m^3
    modelConstants["bio_mx"] = p[1]  # kg/m^3

    modelConstants["K_e"] = p[2]  # mol/m^3
    # modelConstants['K_s'] = p[3] #mol/m^3

    modelConstants["Y_xs"] = p[3]  # g/molS
    modelConstants["Y_ms"] = p[4]  # molB/molS
    modelConstants["Y_os"] = p[5]  # molO/molS

    modelConstants["o2sat"] = 0.214
    modelConstants["K_s"] = 31  # mol/m^3

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


def integratesoln(p, f0, kLa, tfinal, exptdata, ax, nrows, ncols):

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
        row = int((i - 1) / ncols)
        col = int((i - 1) % ncols)
        ax[row][col].set_title(varnames[i])
        ax[row][col].plot(t, f[:, i])
        ax[row][col].plot(exptdata[:, 2 * i - 2], exptdata[:, 2 * i - 1], "ro")


if __name__ == "__main__":

    f0 = np.array([0.214, 0.5, 441.458, 0])
    kLa = 30
    tfinal = 80

    exptdata = np.loadtxt("exptdata_bdo.csv")
    p0 = np.array([17, 11.0, 0.0214, 0.009, 0.88, 0.0467])

    nrows = 2
    ncols = 3
    (fig, ax) = plt.subplots(nrows, ncols)
    integratesoln(p0, f0, kLa, tfinal, exptdata, ax, nrows, ncols)
    plt.show()
