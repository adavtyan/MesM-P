import sys
import random as rnd
from data_stuctures import *
from progressbar import ProgressBar

def place_bead(xyz, xr, nn, box, L, prd):
  for i in range(3):
    xyz[i] += xr[i]
    if prd[i]:
      if xyz[i]<box[i][0]:
        xyz[i] += L[i]
        nn[i] -= 1
      elif xyz[i]>box[i][1]:
        xyz[i] -= L[i]
        nn[i] += 1

  return xyz

def is_outside_domain(xyz, xyz_max, prd):
  if not prd[0] and (xyz[0]<xyz_max[0][0] or xyz[0]>xyz_max[0][1]): return True
  if not prd[1] and (xyz[1]<xyz_max[1][0] or xyz[1]>xyz_max[1][1]): return True
  if not prd[2] and (xyz[2]<xyz_max[2][0] or xyz[2]>xyz_max[2][1]): return True

def add_protein_to_data(c_xyz, prot, data, b_add_to_map=True):
  data.n_mol += 1
  for j in range(prot.N_atom):
    atm = prot.atoms[j]
    xn = c_xyz[j]
    ind = data.n_atoms + atm.id
    mol = data.n_mol
    ty = data.n_atom_types + atm.ty
    data.atoms.append(Atom(ind, mol, ty, xn[0], xn[1], xn[2], xn[3], xn[4], xn[5]))

    # add the atom to the index map
    if b_add_to_map: data.add_coord_to_map(ind)

  for j in range(prot.N_bond):
    bnd = prot.bonds[j]
    ib = data.n_bonds + bnd.id
    ty = data.n_bond_types + bnd.ty
    at1 = data.n_atoms + bnd.at1
    at2 = data.n_atoms + bnd.at2
    data.bonds.append(Bond(ib, ty, at1, at2))

  for j in range(prot.N_ang):
    ang = prot.angles[j]
    ib = data.n_angles + ang.id
    ty = data.n_angle_types + ang.ty
    at1 = data.n_atoms + ang.at1
    at2 = data.n_atoms + ang.at2
    at3 = data.n_atoms + ang.at3
    data.angles.append(Angle(ib, ty, at1, at2, at3))

  for j in range(prot.N_dih):
    dih = prot.dihedrals[j]
    ib = data.n_dihedrals + dih.id
    ty = data.n_dihedral_types + dih.ty
    at1 = data.n_atoms + dih.at1
    at2 = data.n_atoms + dih.at2
    at3 = data.n_atoms + dih.at3
    at4 = data.n_atoms + dih.at4
    data.dihedrals.append(Dihedral(ib, ty, at1, at2, at3))

  data.n_atoms += prot.N_atom
  data.n_bonds += prot.N_bond
  data.n_angles += prot.N_ang
  data.n_dihedrals += prot.N_dih

def add_protein_types(prot, data):
  # Adding bond, angle, and dihedral coefficents
  for ibc in prot.bond_pars:
    data.bond_pars.append(Bond_Par(ibc.type, ibc.kr, ibc.r0, prot.bond_func))

  for iac in prot.angle_pars:
    data.angle_pars.append(Angle_Par(iac.type, iac.kth, iac.th0, prot.angle_func))

  for idc in prot.dih_pars:
    data.dihedral_pars.append(Dihedral_Par(idc.type, idc.kd, idc.nd, idc.dl, prot.dih_func))

  # Incrementing atom, bond, angle, and dihedral types
  data.n_atom_types += prot.N_atom
  data.n_bond_types += prot.N_bond
  data.n_angle_types += prot.N_ang
  data.n_dihedral_types += prot.N_dih


def add_protein_at(prot, data, x0, periodic=[]):
  if periodic==[]: periodic = [True, True, True]

  c_xyz = []
  for j in range(prot.N_atom):
    xb = prot.atoms[j].x
    nn = [0, 0, 0]
    xyz = [x0[0], x0[1], x0[2]]
    place_bead(xyz, xb, nn, data.box, data.L, periodic)

    c_xyz.append([xyz[0], xyz[1], xyz[2], nn[0], nn[1], nn[2]])

  add_protein_to_data(c_xyz, prot, data)

def add_proteins(prot, Nf, data, xyz_max=[], periodic=[], multi=100):
  Np = prot.N_atom
  dmin = prot.dmin
  if xyz_max==[]: xyz_max = data.box
  if periodic==[]: periodic = [True, True, True]
  ext_flags = [1, 1, 1]

  n_tries = 0
  n_out_dom = 0
  n_cur = 0
  # Try up to 100 times more time when the number of filements required
  fr=Nf/100
  if fr<=0: fr=1
  prb = ProgressBar('green', width=50, block='*', empty='o')
  for i in range(Nf*multi):
    n_tries += 1
    outside_domain = False
    
    # Coordinates of center of geometry
    x0 = rnd.uniform(xyz_max[0][0], xyz_max[0][1])
    y0 = rnd.uniform(xyz_max[1][0], xyz_max[1][1])
    z0 = 0.5*(xyz_max[2][0] + xyz_max[2][1])

    # Assigning coordinates for the protein
    c_xyz = []
    for j in range(Np):
      xb = prot.atoms[j].x
      nn = [0, 0, 0]
      xyz = [x0, y0, z0]
      place_bead(xyz, xb, nn, data.box, data.L, periodic)

      # Check if falls outside the domian
      if is_outside_domain(xyz, xyz_max, periodic):
        outside_domain = True
        n_out_dom += 1
        break

      c_xyz.append([xyz[0], xyz[1], xyz[2], nn[0], nn[1], nn[2]])

    if outside_domain or data.has_conflict(c_xyz, dmin, ext_flags): continue

    # if no conflicts add the filament to the data
    n_cur += 1
    add_protein_to_data(c_xyz, prot, data)

    if n_cur%fr==0:
      p = round(100*float(n_cur)/Nf,2)
      prb.render(int(p), 'Adding complex proteins\nNumber of tries %d' % n_tries)

    if n_cur==Nf: break

  # Adding bond, angle, and dihedral coefficents
  for ibc in prot.bond_pars:
    data.bond_pars.append(Bond_Par(ibc.type, ibc.kr, ibc.r0, prot.bond_func))

  for iac in prot.angle_pars:
    data.angle_pars.append(Angle_Par(iac.type, iac.kth, iac.th0, prot.angle_func))

  for idc in prot.dih_pars:
    data.dihedral_pars.append(Dihedral_Par(idc.type, idc.kd, idc.nd, idc.dl, prot.dih_func))
  
  # Incrementing atom, bond, angle, and dihedral types
  data.n_atom_types += Np
  data.n_bond_types += prot.N_bond
  data.n_angle_types += prot.N_ang
  data.n_dihedral_types += prot.N_dih

  p = round(100*float(n_cur)/Nf,2)
  prb.render(int(p), 'Adding complex proteins\nNumber of tries %d' % n_tries)
  print "%d proteins were added out of %d required, with %d attempts made" % (n_cur, Nf, n_tries)
  print "The proteins was generated outside domian %d times\n" % n_out_dom
