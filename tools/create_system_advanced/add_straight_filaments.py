# Author: Aram Davtyan

import sys
import random as rnd
from math import pi, sin, cos, sqrt, ceil, floor
from data_stuctures import *
from progressbar import ProgressBar

C_PI = pi
C_2PI = 2.0*pi
C_PIH = 0.5*pi

THETA_SCAFFOLD_DEFAULT = pi/3.0

THEAT_MIN_DEFAULT = -0.5*pi
THEAT_MAX_DEFAULT = 0.5*pi

ALIGN_NONE = 0
ALIGN_X  = 1
ALIGN_Y  = 2
ALIGN_Z  = 3
ALIGN_XY = 4
ALIGN_XZ = 5
ALIGN_YZ = 6
ALIGN_LAST = 7

TWOD_FLAG_NONE = 0
TWOD_FLAG_XY = 1
TWOD_FLAG_XZ = 2
TWOD_FLAG_YZ = 3
TWOD_FLAG_SP = 4 # On a sphere
TWOD_FLAG_CL = 5 # On a cylinder
TWOD_FLAG_LAST = 6

def extend_xyz(xyz, dx, nn, box, L, prd):
  xyz[0] += dx[0]
  xyz[1] += dx[1]
  xyz[2] += dx[2]
  for i in range(3):
    if prd[i]:
      if xyz[i]<box[i][0]:
        xyz[i] += L[i]
        nn[i] -= 1
      elif xyz[i]>box[i][1]:
        xyz[i] -= L[i]
        nn[i] += 1

  return xyz

def extend_xy(xyz, dx, nn, box, L, prd):
  xyz[0] += dx[0]
  xyz[1] += dx[1]
  for i in range(2):
    if prd[i]:
      if xyz[i]<box[i][0]:
        xyz[i] += L[i]
        nn[i] -= 1
      elif xyz[i]>box[i][1]:
        xyz[i] -= L[i]
        nn[i] += 1

  return xyz

def extend_xi(xyz, dx, nn, box, L, prd, i):
  xyz[i] += dx[i]
  if prd[i]:
    if xyz[i]<box[i][0]:
      xyz[i] += L[i]
      nn[i] -= 1
    elif xyz[i]>box[i][1]:
      xyz[i] -= L[i]
      nn[i] += 1

  return xyz

def extend_xyz_f(xyz, dx, nn, box, L, prd, flag):
  for i in range(3):
    if flag[i]==0: continue
    xyz[i] += dx[i]
    if prd[i]:
      if xyz[i]<box[i][0]:
        xyz[i] += L[i]
        nn[i] -= 1
      elif xyz[i]>box[i][1]:
        xyz[i] -= L[i]
        nn[i] += 1

  return xyz

def add_filament_to_data(c_xyz, fil, data, b_add_to_map=True):
  data.n_mol += 1
  for j in range(fil.Nb):
    data.n_atoms += 1
    xn = c_xyz[j]
    ind = data.n_atoms
    mol = data.n_mol
    #ty = fil.get_type_index(j) + last_atom_type
    ty = fil.get_type_index(j) + data.n_atom_types
    data.atoms.append(Atom(ind, mol, ty, xn[0], xn[1], xn[2], xn[3], xn[4], xn[5]))

    # add the atom to the index map
    if b_add_to_map: data.add_coord_to_map(ind)

    # add bonds, angles, and dihedrals
    if fil.b_bonds and j>0:
      data.n_bonds += 1
      ib = data.n_bonds
      ty = data.n_bond_types
      data.bonds.append(Bond(ib, ty, ind-1, ind))

    if fil.b_angles and j>1:
      data.n_angles += 1
      ib = data.n_angles
      ty = data.n_angle_types
      data.angles.append(Angle(ib, ty, ind-2, ind-1, ind))

    if fil.b_dihedrals and j>2:
      data.n_dihedrals += 1
      ib = data.n_dihedrals
      ty = data.n_dihedral_types
      data.dihedrals.append(Dihedral(ib, ty, ind-3, ind-2, ind-1, ind))

def add_filament_at(fil, data, x0, dx0, periodic=[]):
  # Variable initiation
  Nb = fil.Nb
  R0 = fil.Rb
  L = Nb*R0
  if periodic==[]: periodic = [True, True, True]
  ext_flags = [1, 1, 1]
  ierror = False

  # Make sure that dx0 is unit vector
  norm = dx0[0]*dx0[0] + dx0[1]*dx0[1] + dx0[2]*dx0[2]
  if norm>0.0:
    if norm!=1.0:
      norm_inv = 1.0/norm
      dx0 = [dx[0]*norm_inv, dx[1]*norm_inv, dx[2]*norm_inv]
  else: ierror = True

  if ierror:
    print "add_filament_at() function input error!\n\n"
    sys.exit()

  dx = [R0*dx0[0], R0*dx0[1], R0*dx0[2]]

  c_xyz = []
  xyz = [x0[0], x0[1], x0[2]] 
  nn = [0, 0, 0]
  for j in range(Nb):
    extend_xyz_f(xyz, dx, nn, data.box, data.L, periodic, ext_flags)

    c_xyz.append([xyz[0], xyz[1], xyz[2], nn[0], nn[1], nn[2]])

  add_filament_to_data(c_xyz, fil, data)
    
def gen_helical_scaffold(fil, Nf, data, top_data=[], periodic=True, xyz_max=[], sp_min=0.0):
  # Variable initiation
  Nb = fil.Nb
  R0 = fil.Rb
  Lp = Nb*R0
  dmin = fil.dmin
  #last_atom_type = data.n_atom_types
  #data.n_atom_types += fil.get_ntypes()
  if fil.b_bonds: data.n_bond_types += 1
  if fil.b_angles: data.n_angle_types += 1
  if fil.b_dihedrals: data.n_dihedral_types += 1
  if xyz_max==[]: xyz_max = data.box

  # Sanity checks
  ierror = False
  if len(top_data)!=3: ierror = True
  if ierror:
    print "gen_helical_scaffold() function input error!\n\n"
    sys.exit()

  # Cylinder axis and radius
  Xs = top_data[0]
  Ys = top_data[1]
  Rs = top_data[2]

  # Pre-computing placement parameters
  L_tot = (Lp + sp_min)*Nf # Total length required for Nf proteins with sp spacing between them
  Lz = xyz_max[2][1] - xyz_max[2][0]
  # Length of the helix cannot be shorter than Lz
  if L_tot<=Lz: L_tot = Lz + Lp
  # If the system is periodic, the helix must close on itself accros Z boundary
  if periodic:
    n_turns = ceil(sqrt(L_tot**2 - Lz**2)/(C_2PI*Rs))
    L_tot = sqrt(Lz**2 + (C_2PI*Rs*n_turns)**2)
  # Find the parametric step size dt and vertical step zc
  LL = sqrt(L_tot**2 - Lz**2)
  dt = LL/(Rs*float(Nf))
  zc = Lz*Rs/LL
  # Pre-compute constant portions of dx
  # dx = dx0.{sin(t), cos(t), 1}
  dx0 = [-R0*LL/L_tot, R0*LL/L_tot, Lz/L_tot]

  fr=Nf/100
  if fr<=0: fr=1
  prb = ProgressBar('green', width=50, block='*', empty='o')
  # Generate Nf filaments
  for i in range(Nf):
    # Find a position on the surface of the cylinder
    # This will be used as a midpoint of the filament
    ti = dt*i
    xm = Xs + Rs*cos(ti)
    ym = Ys + Rs*sin(ti)
    zm = xyz_max[2][0] + zc*ti

    # Find displacment and starting position
    dx = [dx0[0]*sin(ti), dx0[1]*cos(ti), dx0[2]]
    x0 = xm - 0.5*Nb*dx[0]
    y0 = ym - 0.5*Nb*dx[1]
    z0 = zm - 0.5*Nb*dx[2]

    # Generate filament bead coordinates
    c_xyz = []
    xyz = [x0, y0, z0]
    nn = [0, 0, 0]
    for j in range(Nb):
      extend_xyz_f(xyz, dx, nn, data.box, data.L, [False, False, periodic], [1, 1, 1])
      c_xyz.append([xyz[0], xyz[1], xyz[2], nn[0], nn[1], nn[2]])

    # adding the filament to the data
    add_filament_to_data(c_xyz, fil, data)

    if (i+1)%fr==0:
      p = round(100.0*float(i+1)/Nf,2)
      prb.render(int(p), 'Adding streight filaments in helical scaffold!')

  # Incrementing atom types
  data.n_atom_types += fil.get_ntypes()

  p = round(100.0,2)
  prb.render(int(p), 'Adding streight filaments in helical scaffold!')
  print "%f filaments were added" % Nf


def gen_circular_scaffold(fil, Nf, data, top_data=[], periodic=True, xyz_max=[], sp_min=0.0, theta=THETA_SCAFFOLD_DEFAULT, complete_turn=False):
  # Variable initiation
  Nb = fil.Nb
  R0 = fil.Rb
  Lp = Nb*R0
  dmin = fil.dmin
  #last_atom_type = data.n_atom_types
  #data.n_atom_types += fil.get_ntypes()
  if fil.b_bonds: data.n_bond_types += 1
  if fil.b_angles: data.n_angle_types += 1
  if fil.b_dihedrals: data.n_dihedral_types += 1
  if xyz_max==[]: xyz_max = data.box

  # Sanity checks
  ierror = False
  if len(top_data)!=3: ierror = True
  if theta==0.0: ierror = True
  if ierror:
    print "gen_circular_scaffold() function input error!\n\n"
    sys.exit()

  # Cylinder axis and radius
  Xs = top_data[0]
  Ys = top_data[1]
  Rs = top_data[2]

  # Pre-compute sin(theta), cos(theta) and placement parameters
  stheta = sin(theta)
  ctheta = cos(theta)
  Rc = Rs/stheta # The longest radius of theta cross-section
  P = C_PI*(3.0*(Rc + Rs) - sqrt((3.0*Rc + Rs)*(Rc + 3.0*Rs))) # Approximate perimeter of the cross-section 
  np = floor(C_2PI*Rc/(Lp + sp_min)) # Number of proteins that fit in one turn
  nt = ceil(float(Nf)/float(np)) # Number of turn required to fit all proteins
  dpsi = C_2PI/np # Angular spacing
  Lz = xyz_max[2][1] - xyz_max[2][0]
  zt = Lz/float(nt) # Spacing between turns

  if complete_turn: Nf = int(np)*int(nt)

  fr=Nf/100
  if fr<=0: fr=1
  prb = ProgressBar('green', width=50, block='*', empty='o')
  # Generate Nf filaments
  for i in range(Nf):
    # Find a position on the surface of the cylinder
    # This will be used as a midpoint of the filament
    iturn = floor(i/np)
    iprot = i%np
    psi = dpsi*iprot
    xm = Xs + Rs*cos(psi)
    ym = Ys + Rs*sin(psi)
    zm = xyz_max[2][0] + zt*iturn - Rc*ctheta*cos(psi)

    # Find displacment and starting position
    z_num = ctheta*sin(psi)/stheta
    norm = 1.0/sqrt(1.0 + z_num**2)
    dx = [-R0*norm*sin(psi), R0*norm*cos(psi), R0*norm*z_num]
    #dx = [R0*stheta*sin(psi), -R0*stheta*cos(psi), R0*ctheta]
    x0 = xm - 0.5*Nb*dx[0]
    y0 = ym - 0.5*Nb*dx[1]
    z0 = zm - 0.5*Nb*dx[2]

    # Generate filament bead coordinates
    c_xyz = []
    xyz = [x0, y0, z0]
    nn = [0, 0, 0]
    for j in range(Nb):
      extend_xyz_f(xyz, dx, nn, data.box, data.L, [False, False, periodic], [1, 1, 1])
      c_xyz.append([xyz[0], xyz[1], xyz[2], nn[0], nn[1], nn[2]])

    # adding the filament to the data
    add_filament_to_data(c_xyz, fil, data)

    if (i+1)%fr==0:
      p = round(100.0*float(i+1)/Nf,2)
      prb.render(int(p), 'Adding streight filaments in circular scaffold!')

  # Incrementing atom types
  data.n_atom_types += fil.get_ntypes()

  p = round(100.0,2)
  prb.render(int(p), 'Adding streight filaments in circular scaffold!')
  print "%f filaments were added" % Nf


def gen_straight_filaments(fil, Nf, data, xyz_max=[], periodic=[], twoD=TWOD_FLAG_NONE, align=ALIGN_NONE, top_data=[], theta_min=THEAT_MIN_DEFAULT, theta_max=THEAT_MAX_DEFAULT, multi=100):
  # Variable initiation
  Nb = fil.Nb
  R0 = fil.Rb
  L = Nb*R0
  dmin = fil.dmin
  #last_atom_type = data.n_atom_types
  #data.n_atom_types += fil.get_ntypes()
  if fil.b_bonds: data.n_bond_types += 1
  if fil.b_angles: data.n_angle_types += 1
  if fil.b_dihedrals: data.n_dihedral_types += 1
  if xyz_max==[]: xyz_max = data.box
  if periodic==[]: periodic = [True, True, True]

  # Sanity checks
  ierror = False
  if twoD<0 or twoD>=TWOD_FLAG_LAST or align<0 or align>=ALIGN_LAST: ierror = True
  if twoD==TWOD_FLAG_XY and (align==ALIGN_Z or align==ALIGN_XZ or align==ALIGN_YZ): ierror = True
  if twoD==TWOD_FLAG_XZ and (align==ALIGN_Y or align==ALIGN_XY or align==ALIGN_YZ): ierror = True
  if twoD==TWOD_FLAG_YZ and (align==ALIGN_X or align==ALIGN_XZ or align==ALIGN_XZ): ierror = True
  if twoD==TWOD_FLAG_SP and len(top_data)!=4: ierror = True
  if twoD==TWOD_FLAG_CL and len(top_data)!=3: ierror = True
  if ierror:
    print "gen_straight_filaments() function input error!\n\n"
    sys.exit()

  # getting the center and radius of the sphere or 
  # center axis and radius of the cylinder, if necassary
  # cylinder is always assumed to be oriented in Z direction 
  # and spanning the box in full Z dimention
  if twoD==TWOD_FLAG_SP:
    Xs = top_data[0]
    Ys = top_data[1]
    Zs = top_data[2]
    Rs = top_data[3]
  elif twoD==TWOD_FLAG_CL:
    Xs = top_data[0]
    Ys = top_data[1]
    Rs = top_data[2]

  # Set alignment based on 2D confinment used
  if twoD==TWOD_FLAG_XY and align==ALIGN_NONE: align = ALIGN_XY
  if twoD==TWOD_FLAG_XZ and align==ALIGN_NONE: align = ALIGN_XZ
  if twoD==TWOD_FLAG_YZ and align==ALIGN_NONE: align = ALIGN_YZ

  # Set extension flags that tell in which coordinates filaments are extended
  ext_flags = [1, 1, 1]
  if twoD==TWOD_FLAG_SP:
    if align==ALIGN_XY: ext_flags = [1, 1, 0]
    elif align==ALIGN_XZ: ext_flags = [1, 0, 1]
    elif align==ALIGN_YZ: ext_flags = [0, 1, 1]
  elif twoD==TWOD_FLAG_CL:
    if align==ALIGN_XY or align==ALIGN_X or align==ALIGN_Y:
      ext_flags = [1, 1, 0]
    elif align==ALIGN_XZ or align==ALIGN_YZ or align==ALIGN_Z:
      ext_flags = [0, 0, 1]
  else:
    if align==ALIGN_X: ext_flags = [1, 0, 0]
    elif align==ALIGN_Y: ext_flags = [0, 1, 0]
    elif align==ALIGN_Z: ext_flags = [0, 0, 1]
    elif align==ALIGN_XY: ext_flags = [1, 1, 0]
    elif align==ALIGN_XZ: ext_flags = [1, 0, 1]
    elif align==ALIGN_YZ: ext_flags = [0, 1, 1]

  # Initializing x0, y0, z0. Will be overwritten if necassary
  x0 = 0.5*(xyz_max[0][0] + xyz_max[0][1])
  y0 = 0.5*(xyz_max[1][0] + xyz_max[1][1])
  z0 = 0.5*(xyz_max[2][0] + xyz_max[2][1])

  n_tries = 0
  n_out_dom = 0
  n_cur = 0
  # Try up to 100 times more time when the number of filements required
  fr=Nf/100
  if fr<=0: fr=1
  prb = ProgressBar('green', width=50, block='*', empty='o')
  for i in range(Nf*multi):
    n_tries += 1

    if twoD==TWOD_FLAG_SP:
      # Finding a random position on the surface of a unite sphere
      # This will be used to find the midpoint of the filament
      psi = rnd.uniform(0, C_2PI)
      spsi = sin(psi)
      cpsi = cos(psi)
      zn = rnd.uniform(-1.0, 1.0)
      znc = sqrt(1.0 - zn*zn)
      xn = cpsi*znc
      yn = spsi*znc
      if align==ALIGN_NONE:
        # Randomly orient the filament on the sphere
        phi = rnd.uniform(0, C_2PI)
        sphi = sin(phi)
        cphi = cos(phi)

        # tau_x = -zn*sin[phi]*cos[psi] - cos[phi]*sin[psi]
        # tau_y = -zn*sin[phi]*sin[psi] + cos[phi]*cos[psi]
        # tau_z = sqrt(1 - zn^2)*sin[phi]
        dx0 = [-zn*sphi*cpsi - cphi*spsi, -zn*sphi*spsi + cphi*cpsi, znc*sphi]
      elif align==ALIGN_X:
        # tau = {sqrt(1 - xn^2), -xn*yn/sqrt(1 - xn^2), -xn*zn/sqrt(1 - xn^2)}
        if xn!=1.0:
          xnc = sqrt(1.0 - xn*xn)
          nnc_inv = xn/xnc
          dx0 = [xnc, -yn*nnc_inv, -zn*nnc_inv]
        else:
          dx0 = [0.0, 0.0, 1.0]
      elif align==ALIGN_Y:
        # tau = {-yn*xn/sqrt(1 - yn^2), sqrt(1 - yn^2), -yn*zn/sqrt(1 - yn^2)}
        if yn!=1.0:
          ync = sqrt(1.0 - yn*yn)
          nnc_inv = yn/ync
          dx0 = [-xn*nnc_inv, ync, -zn*nnc_inv]
        else:
          dx0 = [0.0, 0.0, 1.0]
      elif align==ALIGN_Z:
        # tau = {-zn*xn/sqrt(1 - zn^2), -zn*yn/sqrt(1 - zn^2), sqrt(1 - zn^2)}
        # For this case the same can be achived by {-zn*cos(psi), -zn*cos(psi), sqrt(1 - zn^2)}
        if znc!=0.0:
          nnc_inv = zn/znc
          dx0 = [-xn*nnc_inv, -yn*nnc_inv, znc] 
        else:
          dx0 = [1.0, 0.0, 0.0]
      elif align==ALIGN_XY:
        # Simplifyed expression for tau = {-yn/sqrt(xn^2 + yn^2), xn/sqrt(xn^2 + yn^2), 0.0}
        dx0 = [-sin(psi), cos(psi), 0.0]
      elif align==ALIGN_XZ:
        # tau = {-zn/sqrt(xn^2 + zn^2), 0.0, xn/sqrt(xn^2 + zn^2)}
        norm = sqrt(xn*xn + zn*zn)
        if norm!=0.0:
          invnorm = 1.0/norm
          dx0 = [-zn*invnorm, 0.0, xn*invnorm]
        else:
          dx0 = [0.0, 0.0, 1.0]
      elif align==ALIGN_YZ:
        # tau = {0.0, zn/sqrt(yn^2 + zn^2), -yn/sqrt(yn^2 + zn^2)}
        norm = sqrt(yn*yn + zn*zn)
        if norm!=0.0:
          invnorm = 1.0/norm
          dx0 = [0.0, zn*invnorm, -yn*invnorm]
        else:
          dx0 = [0.0, 0.0, 1.0]

      # Finding the midpoint of the filament
      xm = Xs + Rs*xn
      ym = Ys + Rs*yn
      zm = Zs + Rs*zn
      # Find the starting position
      dx = [R0*dx0[0], R0*dx0[1], R0*dx0[2]]
      x0 = xm - 0.5*Nb*dx[0]
      y0 = ym - 0.5*Nb*dx[1]
      z0 = zm - 0.5*Nb*dx[2]
    elif twoD==TWOD_FLAG_CL:
      # Finding a random position on the surface of the cylinder
      # This will be used as a midpoint of the filament
      psi = rnd.uniform(0, C_2PI)
      xm = Xs + Rs*cos(psi)
      ym = Ys + Rs*sin(psi)
      zm = rnd.uniform(xyz_max[2][0], xyz_max[2][1])
      if align==ALIGN_NONE:
        phi = rnd.uniform(theta_min, theta_max)
        sphi = sin(phi)
        cphi = cos(phi)
        dx0 = [sphi*sin(psi), -sphi*cos(psi), cphi]
      elif align==ALIGN_XY or align==ALIGN_X or align==ALIGN_Y:
        dx0 = [-sin(psi), cos(psi), 0.0]
      elif align==ALIGN_XZ or align==ALIGN_YZ or align==ALIGN_Z:
        dx0 = [0.0, 0.0, 1.0]

      # Finding the starting position
      dx = [R0*dx0[0], R0*dx0[1], R0*dx0[2]]
      x0 = xm - 0.5*Nb*dx[0]
      y0 = ym - 0.5*Nb*dx[1]
      z0 = zm - 0.5*Nb*dx[2]
    else:
      # If filaments are not placed on a sphere or a cylinder
      if twoD!=TWOD_FLAG_YZ:
        x0 = rnd.uniform(xyz_max[0][0], xyz_max[0][1])
      if twoD!=TWOD_FLAG_XZ:
        y0 = rnd.uniform(xyz_max[1][0], xyz_max[1][1])
      if twoD!=TWOD_FLAG_XY:
        z0 = rnd.uniform(xyz_max[2][0], xyz_max[2][1])
    
      if align==ALIGN_NONE:
        phi = rnd.uniform(0, C_2PI)
        theta = rnd.uniform(theta_min, theta_max)
        dx = [R0*sin(theta)*cos(phi), R0*sin(theta)*sin(phi), R0*cos(theta)]
      elif align==ALIGN_X:
        dx = [R0, 0.0, 0.0]
      elif align==ALIGN_Y:
        dx = [0.0, R0, 0.0]
      elif align==ALIGN_Z:
        dx = [0.0, 0.0, R0]
      elif align==ALIGN_XY:
        phi = rnd.uniform(0, 2.0*pi)
        dx = [R0*cos(phi), R0*sin(phi), 0.0]
      elif align==ALIGN_XZ:
        phi = rnd.uniform(0, 2.0*pi)
        dx = [R0*cos(phi), 0.0, R0*sin(phi)]
      elif align==ALIGN_YZ:
        phi = rnd.uniform(0, 2.0*pi)
        dx = [0.0, R0*cos(phi), R0*sin(phi)]

    # pre-check if anything is outside of the allowed domain
    outside_domain = False
    # if placed on sphere or a cylinder also check for the starting point
    if twoD==TWOD_FLAG_SP or twoD==TWOD_FLAG_CL:
      if not periodic[0]:
        if x0<xyz_max[0][0] or x0>xyz_max[0][1]: outside_domain = True
      if not periodic[1]:
        if y0<xyz_max[1][0] or y0>xyz_max[1][1]: outside_domain = True
      if not periodic[2]:
        if z0<xyz_max[2][0] or z0>xyz_max[2][1]: outside_domain = True
    # check for the ending point defined as:
    # xyz = [x0 + RN*dx[0], y0 + RN*dx[1], z0 + RN*dx[2]]
    if not periodic[0]:
      x = x0 + Nb*dx[0]
      if x<xyz_max[0][0] or x>xyz_max[0][1]: outside_domain = True
    if not periodic[1]:
      y = y0 + Nb*dx[1]
      if y<xyz_max[1][0] or y>xyz_max[1][1]: outside_domain = True
    if not periodic[2]:
      z = z0 + Nb*dx[2]
      if z<xyz_max[2][0] or z>xyz_max[2][1]: outside_domain = True
    if outside_domain:
      n_out_dom += 1
      continue

    # main loop where filament coordinates are generated and checked for conflicts
    c_xyz = []
    xyz = [x0, y0, z0]
    nn = [0, 0, 0]
    for j in range(Nb):
      extend_xyz_f(xyz, dx, nn, data.box, data.L, periodic, ext_flags)

      c_xyz.append([xyz[0], xyz[1], xyz[2], nn[0], nn[1], nn[2]])

    if outside_domain or data.has_conflict(c_xyz, dmin, ext_flags): continue

    # if no conflicts add the filament to the data
    n_cur += 1
    add_filament_to_data(c_xyz, fil, data)

    if n_cur%fr==0:
      p = round(100*float(n_cur)/Nf,2)
      prb.render(int(p), 'Adding streight filaments\nNumber of tries %d' % n_tries)

    if n_cur==Nf: break
  
  # Incrementing atom types
  data.n_atom_types += fil.get_ntypes()

  p = round(100*float(n_cur)/Nf,2)
  prb.render(int(p), 'Adding streight filaments\nNumber of tries %d' % n_tries)
  print "%d filaments were added out of %d required, with %d attempts made" % (n_cur, Nf, n_tries)
  print "The filament was generated outside domain %d times\n" % n_out_dom
