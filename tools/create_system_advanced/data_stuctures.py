# Author: Aram Davtyan

import sys
import re
import numpy as np
from math_lib import *

DATA_TY_ATOM = 0
DATA_TY_MOL = 1
DATA_TY_FULL = 2
DATA_TY_EM2 = 3

PROT_AL_NONE = 0
PROT_AL_XY = 1
PROT_AL_X = 2

class Atom:
  id = 0
  mol = 0
  ty = 0
  q = 0
  x = [0.0, 0.0, 0.0]
  quat = [0.0, 0.0, 0.0, 0.0]
  phi = [0.0, 0.0]
  n = [0, 0, 0]

  def __init__(self, id, mol, ty, x, y, z, nx=0, ny=0, nz=0):      
    self.id = id
    self.mol = mol
    self.ty = ty
    self.x = [x, y, z]
    self.n = [nx, ny, nz]
    self.quat = [0.0, 0.0, 0.0, 0.0]
    self.phi = [0.0, 0.0]

  def to_str_atom(self):
    return "%d %d %e %e %e %d %d %d\n" % (self.id, self.ty, self.x[0], self.x[1], self.x[2], self.n[0], self.n[1], self.n[2])
    
  def to_str_mol(self):
    return "%d %d %d %e %e %e %d %d %d\n" % (self.id, self.mol, self.ty, self.x[0], self.x[1], self.x[2], self.n[0], self.n[1], self.n[2])

  def to_str_full(self):
    return "%d %d %d %e %e %e %e %d %d %d\n" % (self.id, self.mol, self.ty, self.q, self.x[0], self.x[1], self.x[2], self.n[0], self.n[1], self.n[2])

  def to_str_em2(self):
    return "%d %d %d %e %e %e %e %e %e %e %e %e %d %d %d\n" % (self.id, self.mol, self.ty, self.quat[0], self.quat[1], self.quat[2], self.quat[3], self.phi[0], self.phi[1], self.x[0], self.x[1], self.x[2], self.n[0], self.n[1], self.n[2])

class Bond:
  id = 0
  ty = 0
  at1 = 0
  at2 = 0

  def __init__(self, id, ty, at1, at2):
    self.id = id
    self.ty = ty
    self.at1 = at1
    self.at2 = at2

  def to_str(self):
    return "%d %d %d %d\n" % (self.id, self.ty, self.at1, self.at2)

class Angle:
  id = 0
  ty = 0
  at1 = 0
  at2 = 0
  at3 = 0

  def __init__(self, id, ty, at1, at2, at3):
    self.id = id
    self.ty = ty
    self.at1 = at1
    self.at2 = at2
    self.at3 = at3

  def to_str(self):
    return "%d %d %d %d %d\n" % (self.id, self.ty, self.at1, self.at2, self.at3)

class Dihedral:
  id = 0
  ty = 0
  at1 = 0
  at2 = 0
  at3 = 0
  at4 = 0

  def __init__(self, id, ty, at1, at2, at3, at4):
    self.id = id
    self.ty = ty
    self.at1 = at1
    self.at2 = at2
    self.at3 = at3
    self.at4 = at4

  def to_str(self):
    return "%d %d %d %d %d %d\n" % (self.id, self.ty, self.at1, self.at2, self.at3, self.at4)

class Bond_Par:
  type = 0
  kr = 0.0
  r0 = 0.0
  func = ""

  def __init__(self, ty, Kr, R0, fn=""):
    self.type = ty
    self.kr = Kr
    self.r0 = R0
    self.func = fn

  def to_str(self):
    if self.func=="":
      return "%d %f %f\n" % (self.type, self.kr, self.r0)
    else:
      return "%d %s %f %f\n" % (self.type, self.func, self.kr, self.r0)

class Angle_Par:
  type = 0
  kth = 0.0
  th0 = 0.0
  func = ""

  def __init__(self, ty, Kth, Th0, fn=""):
    self.type = ty
    self.kth = Kth
    self.th0 = Th0
    self.func = fn

  def to_str(self):
    if self.func=="":
      return "%d %f %f\n" % (self.type, self.kth, self.th0)
    else:
      return "%d %s %f %f\n" % (self.type, self.func, self.kth, self.th0)

class Dihedral_Par:
  type = 0
  kd = 0.0
  nd = 0.0
  dl = 0.0
  func = ""

  def __init__(self, ty, Kd, Nd, Dl, fn=""):
    self.type = ty
    self.kd = Kd
    self.nd = Nd
    self.dl = Dl
    self.func = fn

  def to_str(self):
    if self.func=="":
      return "%d %f %f %f\n" % (self.type, self.kd, self.nd, self.dl)
    else:
      return "%d %s %f %f %f\n" % (self.type, self.func, self.kd, self.nd, self.dl)

class Data:
  data_file_type = DATA_TY_ATOM
  
  n_atoms = 0
  n_bonds = 0
  n_angles = 0
  n_dihedrals = 0
    
  n_atom_types = 0
  n_bond_types = 0
  n_angle_types = 0
  n_dihedral_types = 0

  n_mol = 0

  box = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
  L = [0.0, 0.0, 0.0]
  L_half = [0.0, 0.0, 0.0]
  
  masses = []
  atoms = []
  bonds = []
  angles = []
  dihedrals = []

  bond_pars = []
  angle_pars = []
  dihedral_pars = []

  smaps = [[], [], []]

  def __init__(self, dty, b):
    self.data_file_type = dty
    self.box[0][0] = b[0][0]
    self.box[0][1] = b[0][1]
    self.box[1][0] = b[1][0]
    self.box[1][1] = b[1][1]
    self.box[2][0] = b[2][0]
    self.box[2][1] = b[2][1]

    self.L[0] = b[0][1] - b[0][0]
    self.L[1] = b[1][1] - b[1][0]
    self.L[2] = b[2][1] - b[2][0]

    self.L_half[0] = 0.5*self.L[0]
    self.L_half[1] = 0.5*self.L[1]
    self.L_half[2] = 0.5*self.L[2]

  def write_data_file(self, out):
    out.write("BD data file\n\n")
    out.write("%d atoms\n" % self.n_atoms)
    out.write("%d bonds\n" % self.n_bonds)
    out.write("%d angles\n" % self.n_angles)
    out.write("%d dihedrals\n" % self.n_dihedrals)
    out.write("0 impropers\n\n")
    
    out.write("%d atom types\n" % self.n_atom_types)
    out.write("%d bond types\n" % self.n_bond_types)
    out.write("%d angle types\n" % self.n_angle_types)
    out.write("%d dihedral types\n" % self.n_dihedral_types)
    out.write("0 improper types\n\n")
    
    out.write("%e %e xlo xhi\n" % (self.box[0][0], self.box[0][1]))
    out.write("%e %e ylo yhi\n" % (self.box[1][0], self.box[1][1]))
    out.write("%e %e zlo zhi\n\n" % (self.box[2][0], self.box[2][1]))
    
    out.write("Masses\n\n")
    for i in range(len(self.masses)):
      im = self.masses[i]
      out.write("%d %e\n" % (i+1, im))
    out.write("\n")
    
    out.write("Atoms\n\n")
    for ia in self.atoms:
      str = ""
      if self.data_file_type==DATA_TY_ATOM: str = ia.to_str_atom()
      elif self.data_file_type==DATA_TY_MOL: str = ia.to_str_mol()
      elif self.data_file_type==DATA_TY_FULL: str = ia.to_str_full()
      elif self.data_file_type==DATA_TY_EM2: str = ia.to_str_em2()
      else:
          print "Error: Unknown Data Type\n"
          sys.exit()
      out.write(str)

    if self.n_bonds>0: out.write("\nBonds\n\n")
    for ib in self.bonds:
       out.write(ib.to_str())

    if self.n_angles>0: out.write("\nAngles\n\n")
    for ia in self.angles:
       out.write(ia.to_str())

    if self.n_dihedrals>0: out.write("\nDihedrals\n\n")
    for idi in self.dihedrals:
       out.write(idi.to_str())

    if len(self.bond_pars)>0: out.write("\nBond Coeffs\n\n")
    for ibc in self.bond_pars:
      out.write(ibc.to_str())

    if len(self.angle_pars)>0: out.write("\nAngle Coeffs\n\n")
    for iac in self.angle_pars:
      out.write(iac.to_str())

    if len(self.dihedral_pars)>0: out.write("\nDihedral Coeffs\n\n")
    for idc in self.dihedral_pars:
      out.write(idc.to_str())

  def reset_smaps(self):
    self.smaps = [[], [], []]

  def get_mapped_coords(self, i, j):
    ind = self.smaps[j][i]
    return [self.atoms[ind].x[0], self.atoms[ind].x[1], self.atoms[ind].x[2]]

  def get_mapped_coord(self, i, j, k):
    return self.atoms[self.smaps[j][i]].x[k]
    
  def search_maps(self, x):
    nn = len(self.smaps[0])
    if nn==0: return [0, 0, 0]

    ii = [-1, -1, -1]

    for j in range(3):
      if x[j]>=self.get_mapped_coord(-1, j, j): ii[j] = nn
      elif x[j]<=self.get_mapped_coord(0, j, j): ii[j] = 0
      else:
        i0 = 0
        il = nn - 1
        while (il>i0+1):
          im = int((i0 + il)/2)
          if x[j]<self.get_mapped_coord(im, j, j): il = im
          else: i0 = im
        ii[j] = il

    return ii

  def check_for_conflicts(self, inds, ix, dmin, dmin_sq, dim_flags=[1,1,1]):
    na = len(self.smaps[0])
    ndim = 3
    for j in range(ndim):
      if dim_flags[j]==0: continue
      i1 = inds[j]
      i2 = inds[j] - 1
      for i in range(1,na):
        i1 -= 1
        i2 += 1
        if i1<0: i1 += na
        if i2>=na: i2 -= na

        dx1 = abs(ix[j] - self.get_mapped_coord(i1, j, j))
        while dx1>self.L_half[j]: dx1 = abs(dx1 - self.L[j])
        dx2 = abs(self.get_mapped_coord(i2, j, j) - ix[j])
        while dx2>self.L_half[j]: dx2 = abs(dx2 - self.L[j])

        if dx1 > dmin and dx2 > dmin: break

        if dx1 < dmin:
          ax = self.get_mapped_coords(i1, j)
          dx = [ix[0] - ax[0], ix[1] - ax[1], ix[2] - ax[2]]
          if abs(dx[0])>self.L_half[j]: dx[0] = self.L[j] - abs(dx[0])
          if abs(dx[1])>self.L_half[j]: dx[1] = self.L[j] - abs(dx[1])
          if abs(dx[2])>self.L_half[j]: dx[2] = self.L[j] - abs(dx[2])
          dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]
          if dsq < dmin_sq: return True
        if dx2 < dmin:
          ax = self.get_mapped_coords(i2, j)
          dx = [ix[0] - ax[0], ix[1] - ax[1], ix[2] - ax[2]]
          if abs(dx[0])>self.L_half[j]: dx[0] = self.L[j] - abs(dx[0])
          if abs(dx[1])>self.L_half[j]: dx[1] = self.L[j] - abs(dx[1])
          if abs(dx[2])>self.L_half[j]: dx[2] = self.L[j] - abs(dx[2])
          dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]
          if dsq < dmin_sq: return True

    return False

  def has_conflict(self, xyz, dmin, dim_flags=[1,1,1]):
    if len(self.smaps[0])!=len(self.smaps[1]) or len(self.smaps[1])!=len(self.smaps[2]):
        print "Error: Index Maps have different lengths!\n"
        sys.exit()

    dmin_sq = dmin*dmin
    for ix in xyz:
      ii = self.search_maps(ix)
      if self.check_for_conflicts(ii, ix, dmin, dmin_sq, dim_flags): return True

  def add_coord_to_map(self, ind):
    ii = self.search_maps(self.atoms[ind-1].x)
    self.smaps[0].insert(ii[0], ind-1)
    self.smaps[1].insert(ii[1], ind-1)
    self.smaps[2].insert(ii[2], ind-1)

class Filament:
  Nb = 0
  Rb = 0.0
  dmin = 0.0
  dmin_sq = 0.0
  b_bonds = False
  b_angles = False
  b_dihedrals = False
  ty_map = []

  def __init__(self, nb, rb, dm, b_bd=False, b_ang=False, b_dih=False):
    self.Nb = nb
    self.Rb = rb
    self.dmin = dm
    self.dmin_sq = dm*dm
    self.b_bonds = b_bd
    self.b_angles = b_ang
    self.b_dihedrals = b_dih

  def set_type_map_from_list(self, tmap):
    if len(tmap)!=Nb:
      print "Error: Wring size of filament type map\n"
      sys.exit()
    self.ty_map = tmap

  def set_type_map_lin(self):
    self.ty_map = []
    for i in range(self.Nb): self.ty_map.append(i)

  def set_type_map_sym(self):
    self.ty_map = []
    for i in range(self.Nb):
        ty = i+1
        if i>self.Nb/2 + self.Nb%2 - 1: ty = self.Nb - i 
        self.ty_map.append(ty)

  def get_ntypes(self):
    if self.ty_map==[]: return 1
    else: return max(self.ty_map)

  def get_type_index(self, ind):
    if self.ty_map==[]: return 1
    else: return self.ty_map[ind]

class Protein:
  N_atom = 0
  N_bond = 0
  N_ang = 0
  N_dih = 0
  atoms = []
  bonds = []
  angles = []
  dihedrals = []
  bond_pars = []
  angle_pars = []
  dih_pars = []
  bond_func = ""
  angle_func = ""
  dih_func = ""
  dmin = 0.0
  dmin_sq = 0.0

  def __init__(self, dm, xyz_file="", par_file="", align_n=None, align_n2=None, align_direction=PROT_AL_XY, align_remove_cog=True):
    self.dmin = dm
    self.dmin_sq = dm*dm

    if xyz_file!="": self.read_from_xyz(xyz_file)
    if par_file!="": self.read_par_file(par_file)
    if align_n!=None: self.align(align_n, align_n2, align_direction, align_remove_cog)

  def add_atom(self, id, ty, x, y, z, nx=0, ny=0, nz=0):
    mol = 1
    self.N_atom += 1
    self.atoms.append(Atom(id, mol, ty, x, y, z, nx, ny, nz))

  def add_bond(self, at1, at2, kr, r0):
    self.N_bond += 1
    ty = self.N_bond
    self.bonds.append(Bond(self.N_bond, ty, at1, at2))
    self.bond_pars.append(Bond_Par(ty, kr, r0))

  def add_angle(self, at1, at2, at3, kth, th0):
    self.N_ang += 1
    ty = self.N_ang
    self.angles.append(Angle(self.N_ang, ty, at1, at2, at3))
    self.angle_pars.append(Angle_Par(ty, kth, th0))

  def add_dihedral(self, at1, at2, at3, at4, kd, nd, dl):
    self.N_dih += 1
    ty = self.N_dih
    self.dihedrals.append(Dihedral(self.N_dih, ty, at1, at2, at3, at4))
    self.dih_pars.append(Dihedral_Par(ty, kd, nd, dl))

  def read_from_xyz(self, xyz=""):
    out = open(xyz)
    iline = 0
    nn = 0
    for line in out:
      if iline==0:
        nn = int(line)
      elif iline>1:
        line = line.split()

        if len(line)!=4:
          print "XYZ file format error"
          sys.exit()

        tag = line[0]
        x = float(line[1])
        y = float(line[2])
        z = float(line[3])

        ind = int(re.search(r'\d+', tag).group())

        self.add_atom(ind, ind, x, y, z)
      iline += 1
    out.close()

  def write_xyz(self, xyz):
    out = open(xyz, 'w')
    out.write("%d\n" % self.N_atom)
    out.write("\n")
    for i in range(self.N_atom):
      atm = self.atoms[i]
      out.write("A%d\t%f\t%f\t%f\n" % (atm.id, atm.x[0], atm.x[1], atm.x[2]))
    out.close()

  def read_par_file(self, par=""):
    out = open(par)
    state = ""
    for line in out:
      line = line.strip()

      if line=="" or line[0]=="*" or line[0]=="!" or line[0]=="X": continue

      if line=="BONDS" or line=="ANGLES" or line=="DIHEDRALS" or line=="NONBONDED":
        state=line
      else:
        line = line.split()
        if state=="BONDS":
          if len(line)!=4:
            print "XYZ file format error"
            sys.exit()

          at1 = int(re.search(r'\d+', line[0]).group())
          at2 = int(re.search(r'\d+', line[1]).group())
          kr = float(line[2])
          r0 = float(line[3])

          self.add_bond(at1, at2, kr, r0)
        elif state=="ANGLES":
          if len(line)!=5:
            print "XYZ file format error"
            sys.exit()

          at1 = int(re.search(r'\d+', line[0]).group())
          at2 = int(re.search(r'\d+', line[1]).group())
          at3 = int(re.search(r'\d+', line[2]).group())
          kth = float(line[3])
          th0 = float(line[4])

          self.add_angle(at1, at2, at3, kth, th0)
        elif state=="DIHEDRALS":
          if len(line)!=7:
            print "XYZ file format error"
            sys.exit()

          at1 = int(re.search(r'\d+', line[0]).group())
          at2 = int(re.search(r'\d+', line[1]).group())
          at3 = int(re.search(r'\d+', line[2]).group())
          at4 = int(re.search(r'\d+', line[3]).group())
          kd = float(line[4])
          nd = float(line[5])
          dl = float(line[6])

          self.add_dihedral(at1, at2, at3, at4, kn, nd, dl)
    out.close()

  # Align protein so that n and n2 vectors are parallel to XY plane and/or n vector is parallel to X.
  # Possible values for 'direction' are PROT_AL_XY and PROT_AL_X.
  # In the case of PROT_AL_XY, the plane defined by n and n2 is orienated parallel to XY.
  # In the case of PROT_AL_X, n is additionally oriented parallel to X.
  # If n2 vector is not given only n vector is aligned accordingly.
  def align(self, n, n2=None, direction=PROT_AL_XY, remove_cog=True):
    if direction==PROT_AL_NONE: return
    
    # Senity check
    ierror = False
    if type(n)!=list or len(n)!=3: ierror=True
    if np.linalg.norm(n)==0.0: ierror=True
    if n2!=None and len(n2)!=3: ierror=True
    if n2!=None and (np.linalg.norm(n2)==0.0 or np.dot(n,n2)==0.0): ierror=True
    if direction!=PROT_AL_XY and direction!=PROT_AL_X: ierror=True
    if ierror:
      print "align_along_XY() function input error!\n\n"
      sys.exit()

    # Make sure that n and n2 are normalized vectors
    norm = np.linalg.norm(n)
    if norm!=1.0:
      n[0] /= norm
      n[1] /= norm
      n[2] /= norm
    if n2!=None:
      norm = np.linalg.norm(n2)
      if norm!=1.0:
        n2[0] /= norm
        n2[1] /= norm
        n2[2] /= norm

    if n2==None:
      if direction==PROT_AL_XY:
        # Cosine of rotation angle theta
        cth = np.sign(n[2])*np.sqrt(n[0]*n[0] + n[1]*n[1])
        cth_half = np.sqrt(0.5*(1.0 + cth))
        sth_half = np.sqrt(0.5*(1.0 - cth))
        # Rotation axis
        u = np.array([-n[1], n[0], 0.0]) # Cross[k, n]
        u /= np.linalg.norm(u)

        # Quaternion
        q = [cth_half, sth_half*u[0], sth_half*u[1], sth_half*u[2]]
      else:
        # Cosine of rotation angle theta
        cth = n[0]
        cth_half = np.sqrt(0.5*(1.0 + cth))
        sth_half = np.sqrt(0.5*(1.0 - cth))
        # Rotation axis
        u = np.array([0.0, n[2], -n[1]]) # Cross[n, i]
        u /= np.linalg.norm(u)

        # Quaternion
        q = [cth_half, sth_half*u[0], sth_half*u[1], sth_half*u[2]]
    else:
      # direction->PROT_AL_XY
      # Normal to the plane defined by n and n2
      omega = np.cross(n, n2)
      omega /= np.linalg.norm(omega)

      # Cosine of rotation angle theta
      cth = omega[2] # cos(theta)=omega.k
      cth_half = np.sqrt(0.5*(1.0 + cth))
      sth_half = np.sqrt(0.5*(1.0 - cth))
      # Rotation axis
      u = np.array([omega[1], -omega[0], 0.0]) # Cross[omega, k]
      u /= np.linalg.norm(u)

      # Quaternion
      # q = [cth_half, sth_half*u[0], sth_half*u[1], sth_half*u[2]]
      q = [cth_half, sth_half*u[0], sth_half*u[1], 0.0]

      # If n needs to be alliged with X rotate around Z till aligned
      if direction==PROT_AL_X:
        qmat = quat_to_mat_trans(q)
        nxy = vecmat(n, qmat)
        
        # Cosine of rotation angle theta
        cth = nxy[0] # cos(theta)=nxy.i
        cth_half = np.sqrt(0.5*(1.0 + cth))
        sth_half = np.sqrt(0.5*(1.0 - cth))

        # Quaternion
        qx = [cth_half, 0.0, 0.0, -np.sign(nxy[1])*sth_half]

        # Calculate q = qx**q, assuming that q[3]=0 and qx[1]=qx[2]=0
        q = [q[0]*qx[0], q[1]*qx[0] - q[2]*qx[3], q[2]*qx[0] + q[1]*qx[3], q[0]*qx[3]]

    # Rotation matrix
    qmat = quat_to_mat_trans(q)

    # Calculate the center of geometry
    cog = [0.0, 0.0, 0.0]
    for i in range(self.N_atom):
      atm = self.atoms[i].x
      cog[0] += atm[0]
      cog[1] += atm[1]
      cog[2] += atm[2]
    cog[0] /= self.N_atom
    cog[1] /= self.N_atom
    cog[2] /= self.N_atom

    # Rotate protein around its center of geometry
    xyz = np.zeros(3)
    for i in range(self.N_atom):
      atm = self.atoms[i].x
      xyz[0] = atm[0] - cog[0]
      xyz[1] = atm[1] - cog[1]
      xyz[2] = atm[2] - cog[2]

      xyz = vecmat(xyz, qmat)
      if not remove_cog:
        xyz[0] += cog[0]
        xyz[1] += cog[1]
        xyz[2] += cog[2]

      atm[0] = xyz[0]
      atm[1] = xyz[1]
      atm[2] = xyz[2]
