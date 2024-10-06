import numpy as np

N = 7
X = np.linspace(595353.72, 611949.17, N)
Y = np.linspace(8145346.56, 8163662.48, N)

outF0 = open("MCPL_RxPos" + ".xyz", "w")
for i in range(N):
    for j in range(N):
        outF0.write("{0} {1} {2}".format(str(round(X[i])), str(round(Y[j])), str(0.0)))
        outF0.write("\n")
outF0.close()

print(np.sin(np.deg2rad(30)))
