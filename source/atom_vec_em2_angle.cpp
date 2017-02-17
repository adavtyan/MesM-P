/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aram Davtyan
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "atom_vec_em2_angle.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecEM2Angle::AtomVecEM2Angle(LAMMPS *lmp) : AtomVecEM2(lmp)
{
  molecular = 1;
  bonds_allow = angles_allow = 1;
  forceclearflag = 1;

  comm_x_only = comm_f_only = 0;
  size_forward = 11; // 3 coord + 4 quat + 2 phi + 2 phi_half
  size_reverse = 8; // 3 force + 3 torque + 2 dphi
  size_border = 15;  // 3 coord + 4 tag mol type mask + 4 quat + 2 phi + 2 phi_half
  size_velocity = 6; // 3 vel + 3 angmom
  size_data_atom = 12; // tag + mol + type + 4 quat + 2 phi + 3 coords
  size_data_vel = 7;
  xcol_data = 10;

//  size_data_bonus = 10; 

  atom->molecule_flag = 1;
  atom->angmom_flag = atom->torque_flag = 1;
  mass_type = 1;

  bonus = NULL;
  phi = phi_half = dphi = NULL;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecEM2Angle::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  molecule = memory->grow(atom->molecule,nmax,"atom:molecule");

  nspecial = memory->grow(atom->nspecial,nmax,3,"atom:nspecial");
  special = memory->grow(atom->special,nmax,atom->maxspecial,"atom:special");

  num_bond = memory->grow(atom->num_bond,nmax,"atom:num_bond");
  bond_type = memory->grow(atom->bond_type,nmax,atom->bond_per_atom,
                           "atom:bond_type");
  bond_atom = memory->grow(atom->bond_atom,nmax,atom->bond_per_atom,
                           "atom:bond_atom");

  num_angle = memory->grow(atom->num_angle,nmax,"atom:num_angle");
  angle_type = memory->grow(atom->angle_type,nmax,atom->angle_per_atom,
                            "atom:angle_type");
  angle_atom1 = memory->grow(atom->angle_atom1,nmax,atom->angle_per_atom,
                             "atom:angle_atom1");
  angle_atom2 = memory->grow(atom->angle_atom2,nmax,atom->angle_per_atom,
                             "atom:angle_atom2");
  angle_atom3 = memory->grow(atom->angle_atom3,nmax,atom->angle_per_atom,
                             "atom:angle_atom3");

  angmom = memory->grow(atom->angmom,nmax,3,"atom:angmom");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");

  phi = memory->grow(phi,nmax,2,"atom:phi");
  phi_half = memory->grow(phi_half,nmax,2,"atom:phi_half");
  dphi = memory->grow(dphi,nmax*comm->nthreads,2,"atom:dphi");

  bonus = (Bonus *) memory->srealloc(bonus,nmax*sizeof(Bonus),
                                     "atom:bonus");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecEM2Angle::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  molecule = atom->molecule;
  nspecial = atom->nspecial; special = atom->special;
  num_bond = atom->num_bond; bond_type = atom->bond_type;
  bond_atom = atom->bond_atom;
  num_angle = atom->num_angle; angle_type = atom->angle_type;
  angle_atom1 = atom->angle_atom1; angle_atom2 = atom->angle_atom2;
  angle_atom3 = atom->angle_atom3;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecEM2Angle::copy(int i, int j, int delflag)
{
  int k;

  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  molecule[j] = molecule[i];

  num_bond[j] = num_bond[i];
  for (k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
  }

  num_angle[j] = num_angle[i];
  for (k = 0; k < num_angle[j]; k++) {
    angle_type[j][k] = angle_type[i][k];
    angle_atom1[j][k] = angle_atom1[i][k];
    angle_atom2[j][k] = angle_atom2[i][k];
    angle_atom3[j][k] = angle_atom3[i][k];
  }

  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];
  for (k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];

  angmom[j][0] = angmom[i][0];
  angmom[j][1] = angmom[i][1];
  angmom[j][2] = angmom[i][2];

  phi[j][0] = phi[i][0];
  phi[j][1] = phi[i][1];

  phi_half[j][0] = phi_half[i][0];
  phi_half[j][1] = phi_half[i][1];

  dphi[j][0] = dphi[i][0];
  dphi[j][1] = dphi[i][1];

  copy_bonus(i,j);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2Angle::pack_border(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  double *quat;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      quat = bonus[j].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = phi[j][0];
      buf[m++] = phi[j][1];
      buf[m++] = phi_half[j][0];
      buf[m++] = phi_half[j][1];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      quat = bonus[j].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = phi[j][0];
      buf[m++] = phi[j][1];
      buf[m++] = phi_half[j][0];
      buf[m++] = phi_half[j][1];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2Angle::pack_border_vel(int n, int *list, double *buf,
                                      int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  double *quat;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
      quat = bonus[j].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = phi[j][0];
      buf[m++] = phi[j][1];
      buf[m++] = phi_half[j][0];
      buf[m++] = phi_half[j][1];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = angmom[j][0];
      buf[m++] = angmom[j][1];
      buf[m++] = angmom[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = ubuf(molecule[j]).d;
        quat = bonus[j].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = phi[j][0];
        buf[m++] = phi[j][1];
        buf[m++] = phi_half[j][0];
        buf[m++] = phi_half[j][1];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = ubuf(molecule[j]).d;
        quat = bonus[j].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = phi[j][0];
        buf[m++] = phi[j][1];
        buf[m++] = phi_half[j][0];
        buf[m++] = phi_half[j][1];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2Angle::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(molecule[j]).d;
    quat = bonus[j].quat;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = phi[j][0];
    buf[m++] = phi[j][1];
    buf[m++] = phi_half[j][0];
    buf[m++] = phi_half[j][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecEM2Angle::unpack_border(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    quat = bonus[i].quat;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    phi[i][0] = buf[m++];
    phi[i][1] = buf[m++];
    phi_half[i][0] = buf[m++];
    phi_half[i][1] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecEM2Angle::unpack_border_vel(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    quat = bonus[i].quat;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    phi[i][0] = buf[m++];
    phi[i][1] = buf[m++];
    phi_half[i][0] = buf[m++];
    phi_half[i][1] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    angmom[i][0] = buf[m++];
    angmom[i][1] = buf[m++];
    angmom[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2Angle::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    molecule[i] = (tagint) ubuf(buf[m++]).i;
    quat = bonus[i].quat;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    phi[i][0] = buf[m++];
    phi[i][1] = buf[m++];
    phi_half[i][0] = buf[m++];
    phi_half[i][1] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecEM2Angle::pack_exchange(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  double *quat = bonus[i].quat;
  buf[m++] = quat[0];
  buf[m++] = quat[1];
  buf[m++] = quat[2];
  buf[m++] = quat[3];
  buf[m++] = phi[i][0];
  buf[m++] = phi[i][1];
  buf[m++] = phi_half[i][0];
  buf[m++] = phi_half[i][1];

  buf[m++] = ubuf(molecule[i]).d;

  buf[m++] = ubuf(num_bond[i]).d;
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(bond_type[i][k]).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }

  buf[m++] = ubuf(num_angle[i]).d;
  for (k = 0; k < num_angle[i]; k++) {
    buf[m++] = ubuf(angle_type[i][k]).d;
    buf[m++] = ubuf(angle_atom1[i][k]).d;
    buf[m++] = ubuf(angle_atom2[i][k]).d;
    buf[m++] = ubuf(angle_atom3[i][k]).d;
  }

  buf[m++] = ubuf(nspecial[i][0]).d;
  buf[m++] = ubuf(nspecial[i][1]).d;
  buf[m++] = ubuf(nspecial[i][2]).d;
  for (k = 0; k < nspecial[i][2]; k++) buf[m++] = ubuf(special[i][k]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2Angle::unpack_exchange(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  double *quat = bonus[nlocal].quat;
  quat[0] = buf[m++];
  quat[1] = buf[m++];
  quat[2] = buf[m++];
  quat[3] = buf[m++];
  phi[nlocal][0] = buf[m++];
  phi[nlocal][1] = buf[m++];
  phi_half[nlocal][0] = buf[m++];
  phi_half[nlocal][1] = buf[m++];

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;

  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_angle[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_angle[nlocal]; k++) {
    angle_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    angle_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  nspecial[nlocal][0] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][1] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][2] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < nspecial[nlocal][2]; k++)
    special[nlocal][k] = (tagint) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecEM2Angle::size_restart()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    n += 25 + 2*num_bond[i] + 4*num_angle[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including bonus data
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecEM2Angle::pack_restart(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  buf[m++] = bonus[i].quat[0];
  buf[m++] = bonus[i].quat[1];
  buf[m++] = bonus[i].quat[2];
  buf[m++] = bonus[i].quat[3];
  buf[m++] = phi[i][0];
  buf[m++] = phi[i][1];
  buf[m++] = phi_half[i][0];
  buf[m++] = phi_half[i][1];

  buf[m++] = ubuf(molecule[i]).d;

  buf[m++] = ubuf(num_bond[i]).d;
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(MAX(bond_type[i][k],-bond_type[i][k])).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }

  buf[m++] = ubuf(num_angle[i]).d;
  for (k = 0; k < num_angle[i]; k++) {
    buf[m++] = ubuf(MAX(angle_type[i][k],-angle_type[i][k])).d;
    buf[m++] = ubuf(angle_atom1[i][k]).d;
    buf[m++] = ubuf(angle_atom2[i][k]).d;
    buf[m++] = ubuf(angle_atom3[i][k]).d;
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecEM2Angle::unpack_restart(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  double *quat = bonus[nlocal].quat;
  quat[0] = buf[m++];
  quat[1] = buf[m++];
  quat[2] = buf[m++];
  quat[3] = buf[m++];
  phi[nlocal][0] = buf[m++];
  phi[nlocal][1] = buf[m++];
  phi_half[nlocal][0] = buf[m++];
  phi_half[nlocal][1] = buf[m++];

  molecule[nlocal] = (tagint) ubuf(buf[m++]).i;

  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  num_angle[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_angle[nlocal]; k++) {
    angle_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    angle_atom1[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom2[nlocal][k] = (tagint) ubuf(buf[m++]).i;
    angle_atom3[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecEM2Angle::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  angmom[nlocal][0] = 0.0;
  angmom[nlocal][1] = 0.0;
  angmom[nlocal][2] = 0.0;

  double *quat = bonus[nlocal].quat;
  quat[0] = 0.0;
  quat[1] = 0.0;
  quat[2] = 0.0;
  quat[3] = 0.0;
  phi[nlocal][0] = 0.0;
  phi[nlocal][1] = 0.0;
  phi_half[nlocal][0] = 0.0;
  phi_half[nlocal][1] = 0.0;
  dphi[nlocal][0] = 0.0;
  dphi[nlocal][1] = 0.0;

  molecule[nlocal] = 0;
  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;
  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEM2Angle::data_atom(double *coord, imageint imagetmp, 
                                 char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  molecule[nlocal] = ATOTAGINT(values[1]);
  type[nlocal] = atoi(values[2]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  double *quat = bonus[nlocal].quat;
  quat[0] = atof(values[3]);
  quat[1] = atof(values[4]);
  quat[2] = atof(values[5]);
  quat[3] = atof(values[6]);
  phi[nlocal][0] = atof(values[7]);
  phi[nlocal][1] = atof(values[8]);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  angmom[nlocal][0] = 0.0;
  angmom[nlocal][1] = 0.0;
  angmom[nlocal][2] = 0.0;

  phi_half[nlocal][0] = phi[nlocal][0];
  phi_half[nlocal][1] = phi[nlocal][1];
  dphi[nlocal][0] = 0.0;
  dphi[nlocal][1] = 0.0;

  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecEM2Angle::data_atom_hybrid(int nlocal, char **values)
{
  molecule[nlocal] = ATOTAGINT(values[0]);

  double *quat = bonus[nlocal].quat;
  quat[0] = atof(values[1]);
  quat[1] = atof(values[2]);
  quat[2] = atof(values[3]);
  quat[3] = atof(values[4]);
  phi[nlocal][0] = atof(values[5]);
  phi[nlocal][1] = atof(values[6]);

  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;

  return 7;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecEM2Angle::pack_data(double **buf)
{
  double *quat;

  int m;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    m = 0;
    buf[i][m++] = ubuf(tag[i]).d;
    buf[i][m++] = ubuf(molecule[i]).d;
    buf[i][m++] = ubuf(type[i]).d;

    quat = bonus[i].quat;
    buf[i][m++] = quat[0];
    buf[i][m++] = quat[1];
    buf[i][m++] = quat[2];
    buf[i][m++] = quat[3];
    buf[i][m++] = phi[i][0];
    buf[i][m++] = phi[i][1];

    buf[i][m++] = x[i][0];
    buf[i][m++] = x[i][1];
    buf[i][m++] = x[i][2];
    buf[i][m++] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][m++] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][m++] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecEM2Angle::pack_data_hybrid(int i, double *buf)
{
  double *quat = bonus[i].quat;

  int m = 0;
  buf[m++] = ubuf(molecule[i]).d;
  buf[m++] = quat[0];
  buf[m++] = quat[1];
  buf[m++] = quat[2];
  buf[m++] = quat[3];
  buf[m++] = phi[i][0]; 
  buf[m++] = phi[i][1]; 

  return m; 
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecEM2Angle::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " " TAGINT_FORMAT
            " %d %-1.16e %-1.16e %-1.16e %-1.16e "
            "%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(tagint) ubuf(buf[i][1]).i,(int) ubuf(buf[i][2]).i,
            buf[i][3],buf[i][4],buf[i][5],buf[i][6],
            buf[i][7],buf[i][8],
            buf[i][9],buf[i][10],buf[i][11],
            (int) ubuf(buf[i][12]).i,(int) ubuf(buf[i][13]).i,
            (int) ubuf(buf[i][14]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecEM2Angle::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," " TAGINT_FORMAT " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e",
          (tagint) ubuf(buf[0]).i,buf[1],buf[2],buf[3],buf[4],buf[5],buf[6]);

  return 7;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecEM2Angle::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("angmom")) bytes += memory->usage(angmom,nmax,3);
  if (atom->memcheck("torque")) bytes += memory->usage(torque,nmax*comm->nthreads,3);
  
  if (atom->memcheck("phi")) bytes += memory->usage(phi,nmax,2);
  if (atom->memcheck("phi_half")) bytes += memory->usage(phi_half,nmax,2);
  if (atom->memcheck("dphi")) bytes += memory->usage(dphi,nmax*comm->nthreads,2);

  if (atom->memcheck("molecule")) bytes += memory->usage(molecule,nmax);
  if (atom->memcheck("nspecial")) bytes += memory->usage(nspecial,nmax,3);
  if (atom->memcheck("special"))
    bytes += memory->usage(special,nmax,atom->maxspecial);

  if (atom->memcheck("num_bond")) bytes += memory->usage(num_bond,nmax);
  if (atom->memcheck("bond_type"))
    bytes += memory->usage(bond_type,nmax,atom->bond_per_atom);
  if (atom->memcheck("bond_atom"))
    bytes += memory->usage(bond_atom,nmax,atom->bond_per_atom);

  if (atom->memcheck("num_angle")) bytes += memory->usage(num_angle,nmax);
  if (atom->memcheck("angle_type"))
    bytes += memory->usage(angle_type,nmax,atom->angle_per_atom);
  if (atom->memcheck("angle_atom1"))
    bytes += memory->usage(angle_atom1,nmax,atom->angle_per_atom);
  if (atom->memcheck("angle_atom2"))
    bytes += memory->usage(angle_atom2,nmax,atom->angle_per_atom);
  if (atom->memcheck("angle_atom3"))
    bytes += memory->usage(angle_atom3,nmax,atom->angle_per_atom);

  bytes += nmax*sizeof(Bonus);

  return bytes;
}
