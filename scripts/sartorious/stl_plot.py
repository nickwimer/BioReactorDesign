import argparse
import sys

import numpy as np

from bird.utilities.stl_plotting import plotSTL, plt, pretty_labels

if __name__ == "__main__":
    
    axes = plotSTL("sparger.stl")
    pretty_labels("x", "y", zlabel="z", fontsize=14, ax=axes)
    plt.show()
