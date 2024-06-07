import matplotlib.pyplot as plt
import numpy as np

varnames = ["[biomass]", "[glucose]", "[muconate]", "OUR"]
ylabels = ["mol/m^3", "mol/m^3", "mol/m^3", "mol/m^3/h"]
NVAR = 4

if __name__ == "__main__":

    data = np.loadtxt("timehist.dat")
    t = data[:,0]
    (fig, ax) = plt.subplots(2, 2, constrained_layout=True)
    for i in range(1, NVAR+1):
        row = int((i - 1) / 2)
        col = int((i - 1) % 2)
        ax[row][col].set_title(varnames[i-1])
        ax[row][col].plot(t, data[:,i])
        ax[row][col].set_xlabel("Time (h)")
        ax[row][col].set_ylabel(ylabels[i-1])
    plt.show()
