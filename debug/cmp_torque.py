import sys
import numpy as np

if len(sys.argv)!=3:
  print sys.argv[0], "lammps_file check_file"
  sys.exit()

lfile = sys.argv[1]
cfile = sys.argv[2]

data1 = np.loadtxt(lfile)
f1 = data1[:,8:11]

data2 = np.loadtxt(cfile)
f2 = data2[:,1:]

df = np.abs(f1-f2)

print np.amax(df)
