import numpy as np
import math as m
import scipy.integrate as integrate
from scipy.integrate import quad
from numpy import *
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))
import matplotlib.pyplot as plt
from matplotlib import ticker



def readfiles(f) ->object:
    with open(f, "r") as file:
        for line in file:
            if line != "\n":
                line_count+= 1

            arr = line.split()
