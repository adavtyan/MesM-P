import sys
import random as rn
import numpy as np
from add_solvent_tools import *

# Number of membrane and protein particles
n_mem = 5882
n_prot = 89134
n_tot = n_mem + n_prot

# Is this a vesicular system or not
ves = True

# Generate coordinates and quats or read from old config
gen = False
two_ty_sol = True
conf_file = "config.em2"
Lx = Ly = Lz = 3600.0 # default box dimentions
Lx = Ly = 5000.0
Lz = 4000.0
z0 = Lz/2.0 # Z coordinate of the sheet, for plain membrane
mrad = 76.78
mradsq = mrad*mrad;
srad = 78.8927
sqr3half = 0.866025 # sqrt(3)/2

# data file name
data_file = "data.em2"

R = 0
com = np.array([0.0, 0.0, 0.0])
xyz = []
quat = []
types = []
if gen:
  if ves:
    pass
  else:
    nx = int(round(Lx/mrad))
    ny = int(round(Ly/(sqr3half*mrad)))
    n_mem = nx*ny;

    print nx, ny, n_mem

    for i in range(ny):
      for j in range(nx):
        x0 = 0.0
        y0 = 0.0
        if i%2==1: x0 += 0.5*mrad

        x = x0 + j*mrad
        y = y0 + i*sqr3half*mrad;
        xyz.append([x,y,z0])
        types.append(1)

#        print x, y, xyz[-1]

        quat.append([1.0, 0.0, 0.0, 0.0])

    n_prot = int(round((Lx/srad)*(Ly/srad)*(Lz/srad) - float(n_mem)))

    n_tot = n_mem + n_prot

    print n_mem
    print n_prot

    # Randomly distribute solvent particles relatove to membrane
    box = [[0, Lx], [0, Ly], [0, Lz-mrad]]
    sol_xyz = add_solvent(n_prot, box, 70.0, 3)

    # Add solvent to the membrane system
    # Shift each solvent particle based on membrane position
    for ix in sol_xyz:
      z = ix[2] + z0 + 0.5*mrad
      ty = 2
      if z>Lz:
        z -= Lz
        if two_ty_sol: ty = 3
      xyz.append([ix[0], ix[1], z])
      quat.append([1.0, 0.0, 0.0, 0.0])
      types.append(ty)

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
  types = data[1:n_tot+1,3]

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
    if i>=n_mem and np.linalg.norm(xyz[i] - com)>R:
      phi_b.append(rn.uniform(-1.0, 1.0))
    else:
      if not gen and two_ty_sol and i>=n_mem: types[i] = 3
      phi_b.append(1.0)
else:
  for i in range(n_tot):
    # Membrane composition
    if i<n_mem:
      phi_m.append(rn.uniform(-1.0, 1.0))
    else:
      phi_m.append(0.0)
    # Protein composition
    if i>=n_mem and xyz[i][2]>z0:
      phi_b.append(rn.uniform(-1.0, 1.0))
    else:
      if not gen and two_ty_sol and i>=n_mem: types[i] = 3
      phi_b.append(1.0)

# Set phi averages to zero
sum_phi_m = 0.0
sum_phi_b = 0.0
n_phi_m = 0.0
n_phi_b = 0.0
for i in range(n_tot):
  if i<n_mem:
    sum_phi_m += phi_m[i]
    n_phi_m += 1.0
  if i>=n_mem and ( (ves and np.linalg.norm(xyz[i] - com)>R) or (not ves and xyz[i][2]>z0) ):
    sum_phi_b += phi_b[i]
    n_phi_b += 1.0
sum_phi_m /= n_phi_m
sum_phi_b /= n_phi_b

for i in range(n_tot):
  if i<n_mem: phi_m[i] -= sum_phi_m
  if i>=n_mem and ( (ves and np.linalg.norm(xyz[i] - com)>R) or (not ves and xyz[i][2]>z0) ):
    phi_b[i] -= sum_phi_b

sum_phi_m = 0.0
sum_phi_b = 0.0
for i in range(n_tot):
  sum_phi_m += phi_m[i]
  if i>=n_mem and ( (ves and np.linalg.norm(xyz[i] - com)>R) or (not ves and xyz[i][2]>z0) ):
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

if two_ty_sol: out.write("3 atom types\n")
else: out.write("2 atom types\n")
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
out.write("2 1.0\n")
if two_ty_sol: out.write("3 1.0\n\n")
else: out.write("\n")

# Atoms

# For debuging
n_tot = len(xyz)

out.write("Atoms\n\n")
for i in range(n_tot):
#  ty = 1
#  if i>=n_mem: ty = 2
  ty = types[i]
#  I = 6
  q = quat[i]
  x = xyz[i]
  phm = phi_m[i]
  phb = phi_b[i]
  out.write("%d %d %f %f %f %f %.12f %.12f %f %f %f\n" % (i+1, ty, q[0], q[1], q[2], q[3], phm, phb, x[0], x[1], x[2]))
#  out.write("%d %d %f %f %f %f %f %f %f %f %f %f %f %f\n" % (i+1, ty, I, I, I, q[0], q[1], q[2], q[3], phm, phb, x[0], x[1], x[2]))

out.close()
