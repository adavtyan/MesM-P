import numpy as np

data1 = np.loadtxt("lammps_bend_forces")
f1 = data1[:,5:8]

data2 = np.loadtxt("bend_forces")
f2 = data2[:,1:]


df = np.abs(f1-f2)

print np.amax(df)
