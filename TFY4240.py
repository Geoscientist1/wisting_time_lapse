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


# import cmocean



# This function calculates the Fourier Coefficients Cn
def FourierCoefficient(n, L):
    # we define the integrand of Cn as: V0(x)*sin(n*pi*x/L)
    # integrand = lambda x: (1-((x/L)-0.5)**4)*(sin((n*pi*x)/L))  #The first boundary function as defined by equation(5)
    integrand = lambda x: (
    (sin((5 * pi * x) / L) * (sin((n * pi * x) / L))))  # The second boundary function as defined by equation(6)
    cn = (2 / L) * (1 / sinh(n * pi)) * quad(integrand, 0, L)[0]
    # print(cn)
    return cn


# In[102]:


# This function calculates the Potential V(x,y) as a Fourier Sine Series
def potential_series(L, N, nmax):
    potential = np.zeros((N, N))
    x = linspace(0, L, num=N, endpoint=True)
    y = linspace(0, L, num=N, endpoint=True)
    cn = 1
    temp = 0

    for i in range(N):
        for j in range(N):
            for n in range(1, nmax):
                cn = FourierCoefficient(n, L)
                potential[i, j] = cn * sin((n * pi * x[i]) / L) * sinh((n * pi * y[j]) / L) + temp
                temp = potential[i, j]
                if n == nmax - 1:
                    temp = 0
    return potential


# In[103]:


# ===============================================================================================================================

#    Plotting of V0 and comparison of V0 at y=L with the numerical solution
#    Using fillled contour and 3D surface plots to show the results

# ===============================================================================================================================

# Input arguments
L = 3
N = 50

V0 = np.zeros((N, N))
x = linspace(0, L, num=N, endpoint=True)

for i in range(N):
    V0[i] = (sin((5 * pi * x) / L))

fig1 = plt.figure()
X, Y = np.meshgrid(np.linspace(0, L, N), np.linspace(0, L, N))
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$y$', fontsize=18)
plt.contourf(X, Y, V0)
plt.colorbar()
plt.show()
fig1.savefig("v0_c.pdf")

fig2 = plt.figure()
ax = Axes3D(fig2)
surf = ax.plot_surface(X, Y, V0, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$y$', fontsize=18)
fig2.colorbar(surf, shrink=0.5, aspect=5)
fig2.savefig("v0_s.pdf")
print(V0)

# In[104]:


# ===============================================================================================================================

#     Plotting the electric field as the negative gradient of potential
#     Using "pyplot.quiver" to show the results

# ===============================================================================================================================


# Input arguments
L = 3
N = 15
nmax = 15

potential = potential_series(L, N, nmax)
e_field = np.zeros((N, N))

for i in range(N):
    e_field[i, :] = -np.gradient(potential[i, :])

fig = plt.figure()
X, Y = np.meshgrid(np.linspace(0, L, N), np.linspace(0, L, N))
plt.quiver(Y, X, e_field[0, :], e_field[11, :])
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$y$', fontsize=18)
fig.savefig("e_field_.pdf")

# In[105]:


# ===============================================================================================================================

#      Plotting the numerically solved potential with the listed input arguments
#      Using fillled contour and 3D surface plots to show the results

# ===============================================================================================================================
#      Input arguments
L = 3
N = 150
nmax = 15

potential = potential_series(L, N, nmax)
fig = plt.figure()
X, Y = np.meshgrid(np.linspace(0, L, N), np.linspace(0, L, N))
plt.contourf(Y, X, potential)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$y$', fontsize=18)
plt.colorbar()
fig.savefig('pc_.pdf')

fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(Y, X, potential, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$y$', fontsize=18)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
fig.savefig('ps_.pdf')


# In[106]:


# This function investigates the convergence of Fourier Coefficients Cn
def Conv(L, N, nmax):
    V0_appr = np.zeros((N, N))
    x = linspace(0, L, num=N, endpoint=True)
    cn = 1
    temp = 0

    for i in range(N):
        for n in range(1, nmax):
            cn = FourierCoefficient(n, L)
            V0_appr[i] = cn * sinh((n * pi)) * sin((n * pi * x[i]) / L) + temp
            temp = V0_appr[i]
            if n == nmax - 1:
                temp = 0
    return V0_appr


# In[124]:


# ===============================================================================================================================

#       Convergence analysis of the Fourier Coefficients Cn
#       How many terms n do we have to include to get a good enough convergence towards V0(x)

# ===============================================================================================================================

# Input arguments
L = 3
N = 50
nmax = 7

V0 = np.zeros((N, N))
V0_appr = Conv(L, N, nmax)
x = linspace(0, L, num=N, endpoint=True)

for i in range(N):
    V0[i] = (sin((5 * pi * x) / L))

fig = plt.figure()
X, Y = np.meshgrid(np.linspace(0, L, N), np.linspace(0, L, N))
plt.contourf(Y, X, V0_appr)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$y$', fontsize=18)
plt.colorbar()
fig.savefig('convc7_.pdf')

fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(Y, X, V0_appr, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$y$', fontsize=18)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
fig.savefig('convs3_.pdf')

diff = V0 - V0_appr
print(diff)

# In[ ]:




