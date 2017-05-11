# This script allows generation the toplogy file
# for LAMMPS simulations of MesM-P system with explicit I-BARs and actin filaments
# Either flat membrane or vesicular systems can be created

import sys
import random as rn
import numpy as np
from data_stuctures import *
from add_streight_filaments import *
from add_complex_proteins import *
from add_solvent_tools import *

SMALL = 0.000001

MT_VESICLE    = 0
MT_FLAT_SHEET = 1
MT_CYLINDER   = 2

# Main switch keys
mem_topology = MT_VESICLE   # Topology of the system
#mem_topology = MT_CYLINDER   # Topology of the system
#mem_topology = MT_FLAT_SHEET   # Topology of the system
b_solvent = False  # Include solvent or not
two_ty_sol = False   # Use two type of solvents, where one cannot have protein concetration
b_ibar = True # Include I-BAR or not
b_filaments = False # Include actin filaments or not
b_proteins = False # Include complex proteins
data_file = "data.em2"   # data file name
sqr3half = 0.866025 # sqrt(3)/2

# Box dimentions
Lx = Ly = 3600.0 
Lz = 3600.0
L = [Lx, Ly, Lz]
box = [[0, Lx], [0, Ly], [0, Lz]]

# Masses
mem_mass = 1          # mass of the membrane particle
sol_mass = 1          # mass of the solvent particle
ibar_mass = 1         # mass of the ibar particle
actin_mass = 1        # mass of the actin filament particle
prot_mass = 1         # mass of the protein particle

# Solvent parameters
srad = 39.4464      # The size of the solvent particles
srad_min = 35.0     # Minimal distance between two solvent particles
# Standard value for mrad: 76.78
# Standard value for srad_min: 70.0

# Defining parameters for flat and vesicular systems
if mem_topology==MT_VESICLE:
  n_mem = 5882        # Number of mebrane particles
  mrad = 73.9542      # The size of the membrane particles along the surface 
  xM0 = Lx/2.0        # Center of the vesicle
  yM0 = Ly/2.0
  zM0 = Lz/2.0
  # The number of solvent particles will be calculated automatically
  # Standard value for mrad: 73.9542
elif mem_topology==MT_FLAT_SHEET:
  mrad = 35.0         # The size of the membrane particles
  mradsq = mrad*mrad;
  zM0 = Lz/2.0        # z0 coordinte for membrane plane
  # The number of membrane and solvent particles will be calculated automatically
  # Standard value for mrad: 76.78
elif mem_topology==MT_CYLINDER:
  mrad = 35.0         # The size of the membrane particles
  mradsq = mrad*mrad;
  xM0 = Lx/2.0        # x0 coordinte for cylinder axis
  yM0 = Ly/2.0        # y0 coordinte for cylinder axis
  zM0 = Lz/2.0        # Not needed in this case, defined here for generality
  D0 = 1000.0         # Diameter of the cylinder 
  # The number of membrane and solvent particles will be calculated automatically
  # Standard value for mrad: 76.78
else:
  print "Error. Unknown membrane topology\n"
  sys.exit()

# I-BAR parameters
Ni = 500   # Number of I-BARs to add
Nib = 5    # Beads per I-BAR
Rib = 30   # Seperation between I-BAR beads
di_min = 1.5*Rib # Minimum seperation between I-BAR beads
zi0 = zM0 - 2.0*Rib # The plane where I-BARs will be initiated

# Actin filament parameters
Nf = 10    # Number of actin filaments to add
Nfb = 30   # Beads per actin filament
Rfb = 100  # Seperation between actin beads
df_min = 1.5*Rfb # Minimum seperation between actin beads
zfh = zM0 - Rfb/2.0 # Maximum Z coordinates for actin filaments
zfl = 0    # Minimum Z coordinates for actin filaments

# Complex protein parameters
Np = 5               # Number of proteins
dp_min = mrad        # Minimum seperation between proteins
zp0 = zM0 - 2.0*Rib  # The plane where proteins will be initiated
xyz_file = "prot.xyz"
par_file = "prot.par"

# Creating data stuctura
data = Data(DATA_TY_EM2, box)

# Generating membrane coordinates
data.n_atom_types += 1
data.n_mol += 1
data.masses.append(mem_mass)
mem_type = data.n_atom_types
if mem_topology==MT_VESICLE:
  print "Generating a membrane vesicle ..."
  # Finding vesicle radius based on m_mem and diameter of particles 
  R = np.sqrt(n_mem * mrad*mrad / (4.0*np.pi))
  print "R: ", R
  if R>0.5*Lx or R>0.5*Ly or R>0.5*Lz:
    print "Warrning: R is smaller than one of the dimentions!\n"

  # Generate equially spaced vesicle using Fibonacci lattice method
  golden_angle = np.pi * (3.0 - np.sqrt(5.0))
  theta = golden_angle * np.arange(n_mem)
  zu = np.linspace(1.0 - 1.0 / n_mem, 1.0 / n_mem - 1.0, n_mem)
  ru = np.sqrt(1.0 - zu * zu)
  xu = ru * np.cos(theta)
  yu = ru * np.sin(theta)

  for i in range(n_mem):
    # Finding coordinates
    x = R * xu[i] + xM0
    y = R * yu[i] + yM0
    z = R * zu[i] + zM0

    # Addition new atom
    data.n_atoms += 1
    ind = data.n_atoms
    mol = data.n_mol
    ty = data.n_atom_types
    data.atoms.append(Atom(ind, mol, ty, x, y, z))
    
    # Finding quaternion
    z1 = np.sqrt(0.5*(1.0 + zu[i]))
    xy2 = xu[i]**2 + yu[i]**2
    qw = z1
    qx = yu[i]*z1*(zu[i] - 1.0)/xy2
    qy = -xu[i]*z1*(zu[i] - 1.0)/xy2
    qz = 0.0
    data.atoms[-1].quat = [qw, qx, qy, qz]

    # Set membrane and protein composition
    phi_m = rn.uniform(-1.0, 1.0)
    phi_b = 0.0
    data.atoms[-1].phi = [phi_m, phi_b]
elif mem_topology==MT_FLAT_SHEET:
  # First generate a flat sheet membrane
  print "Generating a flat membrane ..."

  # Find the lattice dimentions
  nx = int(round(Lx/mrad))
  ny = int(round(Ly/(sqr3half*mrad)))
  n_mem = nx*ny

  print "nx: %d ny: %d" % (nx, ny)
  print "n_mem:", n_mem
  print

  # Place the membrane particles onto the lattice
  for i in range(ny):
    for j in range(nx):
      x0 = 0.0
      y0 = 0.0
      if i%2==1: x0 += 0.5*mrad

      x = x0 + j*mrad
      y = y0 + i*sqr3half*mrad
      z = zM0

      # Addition new atom
      data.n_atoms += 1
      ind = data.n_atoms
      mol = data.n_mol
      ty = data.n_atom_types
      data.atoms.append(Atom(ind, mol, ty, x, y, z))

      # Set quaternion
      data.atoms[-1].quat = [1.0, 0.0, 0.0, 0.0]

      # Set membrane and protein composition
      phi_m = rn.uniform(-1.0, 1.0)
      phi_b = 0.0
      data.atoms[-1].phi = [phi_m, phi_b]
elif mem_topology==MT_CYLINDER:
  # First generate a cylindrical membrane
  print "Generating a cylindrical membrane ..."

  # Find the lattice dimentions
  nx = int(round(np.pi*D0/mrad))
  ny = int(round(Lz/(sqr3half*mrad)))
  alpha0 = 2.0*np.pi/nx
  R = 0.5*D0
  n_mem = nx*ny

  print "nxy: %d nz: %d" % (nx, ny)
  print "n_mem:", n_mem
  print

  # Place the membrane particles onto the lattice
  for i in range(ny):
    for j in range(nx):
      alpha = j*alpha0
      if i%2==1: alpha += 0.5*alpha0 
      x = R*np.cos(alpha) + xM0
      y = R*np.sin(alpha) + yM0
      z = i*sqr3half*mrad

      # Addition new atom
      data.n_atoms += 1
      ind = data.n_atoms
      mol = data.n_mol
      ty = data.n_atom_types
      data.atoms.append(Atom(ind, mol, ty, x, y, z))

      # Set quaternion
      if fabs(alpha-np.pi)>SMALL:
        qw = 0.5*np.sqrt(1.0 + np.cos(alpha))
        qx = -qw*np.tan(0.5*alpha)
        qy = qw
        qz = -qx
      else:
        qw = qy = 0.0
        qx = 1.0/np.sqrt(2.0)
        qz = -1.0/np.sqrt(2.0)
      data.atoms[-1].quat = [qw, qx, qy, qz]

      # Set membrane and protein composition
      phi_m = rn.uniform(-1.0, 1.0)
      phi_b = 0.0
      data.atoms[-1].phi = [phi_m, phi_b]
else:
  print "Error. Unknown membrane topology\n"
  sys.exit()

# Zero average lipid composition and protein concentration on the membrane
# First calculate the average values
sum_phi_m = 0.0
sum_phi_b = 0.0
n_phi = 0.0
for i in range(data.n_atoms):
  if data.atoms[i].ty==mem_type:
    sum_phi_m += data.atoms[i].phi[0]
    sum_phi_b += data.atoms[i].phi[1]
    n_phi += 1.0
if n_phi!=0:
  sum_phi_m /= n_phi
  sum_phi_b /= n_phi
# Substract the average value
for i in range(data.n_atoms):
  if data.atoms[i].ty==mem_type:
    data.atoms[i].phi[0] -= sum_phi_m
    data.atoms[i].phi[1] -= sum_phi_b
# Calculate the sum again and print them
sum_phi_m = 0.0
sum_phi_b = 0.0
for i in range(data.n_atoms):
  if data.atoms[i].ty==mem_type:
    sum_phi_m += data.atoms[i].phi[0]
    sum_phi_b += data.atoms[i].phi[1]
print "Total membrane phi_m and phi_b:", sum_phi_m, sum_phi_b
print


# Generate I-BAR coordinates
if b_ibar:
  print "Adding I-BAR proteins ..."
  ibar_fil = Filament(Nib, Rib, di_min, True, True)
  ibar_fil.set_type_map_sym()
  
  if mem_topology==MT_VESICLE:
    # Example of placing proteins outside a vesicle
    # The proteins alignmnet parallel to XY plane is used here
    # The following alignments may be used:
    # ALIGN_NONE, ALIGN_X, ALIGN_Y, ALIGN_Z, ALIGN_XY, ALIGN_XZ, ALIGN_YZ
    xyz_max = [box[0], box[1], box[2]]
    gen_straight_filaments(ibar_fil, Ni, data, xyz_max, twoD=TWOD_FLAG_SP, top_data=[xM0, yM0, zM0, R+1.5*Rib], align=ALIGN_XY)
  elif mem_topology==MT_CYLINDER:
    # Example of placing proteins inside a cylinder
    # The cylinder is assumed to be oriented along Z axis
    # The proteins alignmnet parallel to Z axis is used here
    # The following alignments may be used:
    # ALIGN_NONE, ALIGN_Z, ALIGN_XY
    xyz_max = [box[0], box[1], box[2]]
    gen_straight_filaments(ibar_fil, Ni, data, xyz_max, twoD=TWOD_FLAG_CL, top_data=[xM0, yM0, R-1.5*Rib], align=ALIGN_Z)
  elif mem_topology==MT_FLAT_SHEET:
    # Example of placing proteins on a flat membrane
    # The membrane is assumed to be oriented in XY plane
    # The proteins are placed "below" the membrane on Z=zi0 plane
    # and are randonly oriented in XY plane
    # The following alignments may be used for twoD=TWOD_FLAG_XY:
    # ALIGN_NONE, ALIGN_X, ALIGN_Y, ALIGN_XY - equivalent to ALIGN_NONE
    xyz_max = [box[0], box[1], [zi0, zi0]]
    gen_straight_filaments(ibar_fil, Ni, data, xyz_max, twoD=TWOD_FLAG_XY)
  else:
    # Example of placing proteins "below" the membrane
    # The proteins are aligned parallel to XY plane
    # The following alignments may be used:
    # ALIGN_NONE, ALIGN_X, ALIGN_Y, ALIGN_Z, ALIGN_XY, ALIGN_XZ, ALIGN_YZ
    xyz_max = [box[0], box[1], [box[2][0], zi0]]
    gen_straight_filaments(ibar_fil, Ni, data, xyz_max, periodic=[True, True, False], align=ALIGN_XY)

  nty = ibar_fil.get_ntypes()
  for i in range(nty):
    data.masses.append(ibar_mass)

# Generate actin filaments
if b_filaments:
  print "Adding actin filaments ..."
  actin_fil = Filament(Nfb, Rfb, df_min, True, True)
  xyz_max = [box[0], box[1], [zfl, zfh]]
  prd = [True, True, False]
  
  # Reset index maps because filaments are depositied in a different domian from membrane and I-BAR
  data.reset_smaps()
  
  th_max = np.pi/6.0
  gen_straight_filaments(actin_fil, Nf, data, xyz_max, prd, theta_max=th_max)
  data.masses.append(actin_mass)

# Place complex proteins
if b_proteins:
  print "Adding complex proteins ..."
  vec1 = [-10.583, 26.972, 109.369] # Main axis of the protein
  vec2 = [-1.431, 13.432, -8.767] # Secondary axis of the protein
  # Create protein object where vec1 is aligned with X, and plane defined by vec1 and vec2 parelle to XY
  prot = Protein(dp_min, xyz_file, par_file, vec1, vec2, align_direction=PROT_AL_X)

  xyz_max = [box[0], box[1], [zp0, zp0]]
  add_proteins(prot, Np, data, xyz_max)

  for i in range(prot.N_atom):
    data.masses.append(prot_mass)

# Generating solvent coordinates
if b_solvent:
  print "Adding solvent ..."
  # adding type(s) and mass(es) for solvent
  first_sol_type = data.n_atom_types + 1
  if not two_ty_sol:
    data.n_atom_types += 1
    data.masses.append(sol_mass)
  else:
    data.n_atom_types += 2
    data.masses.append(sol_mass)
    data.masses.append(sol_mass)

  n_sol = int(round((Lx/srad)*(Ly/srad)*(Lz/srad) - float(n_mem)))
  print "n_sol:", n_sol

  # Randomly distribute solvent particles in the box
  # For flat sheet reduce the size of the box in Z direction by the radius of one particle
  if mem_topology==MT_FLAT_SHEET:
    box_sol = [[0, Lx], [0, Ly], [0, Lz-mrad]] 
  else:
    box_sol = [[0, Lx], [0, Ly], [0, Lz]]
  sol_xyz = add_solvent(n_sol, box_sol, srad_min, 3)

  # Generate solvent particle type and phi_b
  types = []
  phi_b = []
  sum_phi_b = 0.0
  n_phi_b = 0.0
  for ix in sol_xyz:
    # Shifting Z coordinate in case of flat sheet
    if mem_topology==MT_FLAT_SHEET:
      ix[2] += zM0 + 0.5*mrad
      if ix[2]>Lz: ix[2] -= Lz

    # Assign a correct type to the solvent bead
    ty = first_sol_type
    if two_ty_sol:
      if mem_topology==MT_FLAT_SHEET and ix[2]<zM0:
        ty = first_sol_type + 1
      elif mem_topology==MT_VESICLE and (ix[0]-xM0)**2 + (ix[1]-yM0)**2 + (ix[2]-zM0)**2<R**2:
        ty = first_sol_type + 1
      elif mem_topology==MT_CYLINDER and (ix[0]-xM0)**2 + (ix[1]-yM0)**2<R**2:
        ty = first_sol_type + 1

    # Assigning phi_b and summing over it
    if ty==first_sol_type:
      phi_b_one = rn.uniform(-1.0, 1.0)
      sum_phi_b += phi_b_one
      n_phi_b += 1.0

    types.append(ty)
    phi_b.append(phi_b_one)

  if n_phi_b!=0.0: sum_phi_b /= n_phi_b

  # Add solvent to the membrane system
  for i in range(len(sol_xyz)):
    ix = sol_xyz[i]
    data.n_atoms += 1
    ind = data.n_atoms
    ty = types[i]
    mol = 0
    phi_m_i = 0.0
    phi_b_i = phi_b[i]
    if ty==first_sol_type:  phi_b_i -= sum_phi_b

    # Adding new atom
    data.atoms.append(Atom(ind, mol, ty, ix[0], ix[1], ix[2]))
    data.atoms[-1].quat = [1.0, 0.0, 0.0, 0.0]
    data.atoms[-1].phi = [phi_m_i, phi_b_i]

# write data file

out = open(data_file, 'w')

data.write_data_file(out)

out.close()
