import sys

DATA_TY_ATOM = 0
DATA_TY_MOL = 1
DATA_TY_FULL = 2
DATA_TY_EM2 = 3

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
