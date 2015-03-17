import sys
import random as rn
import numpy as np

# Number of membrane and protein particles
n_mem =5882
n_prot = 89134
n_tot = n_mem + n_prot

# Is this a vesicular system or not
ves = True

# Generate coordinates and quats or read from old config
gen = False
conf_file = "config.em2"
Lx = Ly = Lz = 3600.0 # default box dimentions

# data file name
data_file = "data.em2"

xyz = []
quat = []
if gen:
  pass
else:
  lconv = 10.0 # Distance conversion factor
  data = np.loadtxt(conf_file)
  if data.shape!=(2*n_tot+1, 4):
    print "Config file format error!\n"
    sys.exit()

  # Check to make sure that the types in config 
  # correspond to system definition
  for i in range(n_tot):
    ty = data[i+1,3]
    if (i<n_mem and ty!=1) or (i>=n_mem and ty!=2):
      print "Config file format error!\n"
      sys.exit()

  Lx = lconv*data[0,0]
  Ly = lconv*data[0,1]
  Lz = lconv*data[0,2]

  xyz = lconv*data[1:n_tot+1,:3]
  quat = data[n_tot+1:]

print "Box dimentions" 
print Lx, Ly, Lz
print

# Generate initial phi_m and phi_b compositions
# If ves=True make sure protein composition is absente inside the vesical

phi_m = []
phi_b = []
if ves:
  # Membrane center of mass
  com = np.array([0.0, 0.0, 0.0])
  for i in range(n_mem):
    com += xyz[i,0]
  com /= n_mem

  R = 0.0
  for i in range(n_mem):
    dr = xyz[i] - com
    R += np.linalg.norm(xyz[i] - com)
  R /= n_mem

  for i in range(n_tot):
    # Membrane composition
    if i<n_mem:
      phi_m.append(rn.uniform(-1.0, 1.0))
    else:
      phi_m.append(0.0)
    # Protein composition
    if i<n_mem or np.linalg.norm(xyz[i] - com)>R:
      phi_b.append(rn.uniform(-1.0, 1.0))
    else:
      phi_b.append(1.0)
else:
  for i in range(n_tot):
    # Membrane composition
    if i<n_mem:
      phi_m.append(rn.uniform(-1.0, 1.0))
    else:
      phi_m.append(0.0)
    # Protein composition
    phi_b.append(rn.uniform(-1.0, 1.0))

# Set phi averages to zero
sum_phi_m = 0.0
sum_phi_b = 0.0
for i in range(n_tot):
  if i<n_mem: sum_phi_m += phi_m[i]
  sum_phi_b += phi_b[i]
sum_phi_m /= n_mem
sum_phi_b /= n_tot

for i in range(n_tot):
  if i<n_mem: phi_m[i] -= sum_phi_m
  phi_b[i] -= sum_phi_b

sum_phi_m = 0.0
sum_phi_b = 0.0
for i in range(n_tot):
  sum_phi_m += phi_m[i]
  sum_phi_b += phi_b[i]
print sum_phi_m, sum_phi_b

# write data file

out = open(data_file, 'w')

# Data file header
out.write("EM2 data file\n\n")

out.write("%d atoms\n" % n_tot)
out.write("0 bonds\n")
out.write("0 angles\n")
out.write("0 dihedrals\n")
out.write("0 impropers\n\n")

out.write("2 atom types\n")
out.write("0 bond types\n")
out.write("0 angle types\n")
out.write("0 dihedral types\n")
out.write("0 improper types\n\n")

# Box bounds

out.write("%f %f xlo xhi\n" % (0.0, Lx))
out.write("%f %f ylo yhi\n" % (0.0, Ly))
out.write("%f %f zlo zhi\n\n" % (0.0, Lz))

# Masses

out.write("Masses\n\n")
out.write("1 1.0\n")
out.write("2 1.0\n\n")

# Atoms

out.write("Atoms\n\n")
for i in range(n_tot):
  ty = 1
  if i>=n_mem: ty = 2
#  I = 6
  q = quat[i]
  x = xyz[i]
  phm = phi_m[i]
  phb = phi_b[i]
  out.write("%d %d %f %f %f %f %.12f %.12f %f %f %f\n" % (i+1, ty, q[0], q[1], q[2], q[3], phm, phb, x[0], x[1], x[2]))
#  out.write("%d %d %f %f %f %f %f %f %f %f %f %f %f %f\n" % (i+1, ty, I, I, I, q[0], q[1], q[2], q[3], phm, phb, x[0], x[1], x[2]))

out.close()
