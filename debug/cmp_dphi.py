import sys
import numpy as np

if len(sys.argv)!=3:
  print sys.argv[0], "lammps_file check_file"
  sys.exit()

lfile = sys.argv[1]
cfile = sys.argv[2]

data1 = np.loadtxt(lfile)
f1 = data1[:,11:13]

data2 = np.loadtxt(cfile)
f2 = data2[:,1:]

df = np.abs(f1-f2)

print np.amax(df)

print np.amax(np.abs(f1[:,0]-f2[:,0]))
print np.amax(np.abs(f1[:,1]-f2[:,1]))

#for i in range(len(df)):
#  if df[i][1]>0.0001:
#    print i+1, df[i]
#    print f1[i], f2[i]

