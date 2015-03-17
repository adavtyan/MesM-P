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
#include "atom_vec_em2.h"
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

AtomVecEM2::AtomVecEM2(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  forceclearflag = 1;

  comm_x_only = comm_f_only = 0;
  size_forward = 11; // 3 coord + 4 quat + 2 phi + 2 phi_half
  size_reverse = 8; // 3 force + 3 torque + 2 dphi
  size_border = 14;  // 3 coord + 3 tag type mask + 4 quat + 2 phi + 2 phi_half
  size_velocity = 6; // 3 vel + 3 angmom
  size_data_atom = 11; // tag + type + 4 quat + 2 phi + 3 coords
  size_data_vel = 7;
  xcol_data = 9;

//  size_data_bonus = 10; 

  atom->angmom_flag = atom->torque_flag = 1;
  mass_type = 1;

  bonus = NULL;
  phi = phi_half = dphi = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecEM2::~AtomVecEM2()
{
  memory->sfree(bonus);
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecEM2::grow(int n)
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

void AtomVecEM2::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecEM2::copy(int i, int j, int delflag)
{
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

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset ellipsoid that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecEM2::copy_bonus(int i, int j)
{
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::pack_comm(int n, int *list, double *buf,
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
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
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
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::pack_comm_vel(int n, int *list, double *buf,
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
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
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
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    quat = bonus[j].quat;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = phi[i][0];
    buf[m++] = phi[i][1];
    buf[m++] = phi_half[i][0];
    buf[m++] = phi_half[i][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecEM2::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
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
}

/* ---------------------------------------------------------------------- */

void AtomVecEM2::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
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
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
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

/* ---------------------------------------------------------------------- */

int AtomVecEM2::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
    buf[m++] = dphi[i][0];
    buf[m++] = dphi[i][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
    buf[m++] = dphi[i][0];
    buf[m++] = dphi[i][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecEM2::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
    dphi[j][0] += buf[m++];
    dphi[j][1] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
    dphi[j][0] += buf[m++];
    dphi[j][1] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::pack_border(int n, int *list, double *buf,
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

int AtomVecEM2::pack_border_vel(int n, int *list, double *buf,
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

int AtomVecEM2::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
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

void AtomVecEM2::unpack_border(int n, int first, double *buf)
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

void AtomVecEM2::unpack_border_vel(int n, int first, double *buf)
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

int AtomVecEM2::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
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

int AtomVecEM2::pack_exchange(int i, double *buf)
{
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

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecEM2::unpack_exchange(double *buf)
{
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

int AtomVecEM2::size_restart()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    n += 22;

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

int AtomVecEM2::pack_restart(int i, double *buf)
{
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

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecEM2::unpack_restart(double *buf)
{
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

void AtomVecEM2::create_atom(int itype, double *coord)
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

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEM2::data_atom(double *coord, imageint imagetmp, 
                                 char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  double *quat = bonus[nlocal].quat;
  quat[0] = atof(values[2]);
  quat[1] = atof(values[3]);
  quat[2] = atof(values[4]);
  quat[3] = atof(values[5]);
  phi[nlocal][0] = atof(values[6]);
  phi[nlocal][1] = atof(values[7]);

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

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecEM2::data_atom_hybrid(int nlocal, char **values)
{
  double *quat = bonus[nlocal].quat;
  quat[0] = atof(values[3]);
  quat[1] = atof(values[4]);
  quat[2] = atof(values[5]);
  quat[3] = atof(values[6]);
  phi[nlocal][0] = atof(values[7]);
  phi[nlocal][1] = atof(values[8]);

  return 9;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecEM2::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  angmom[m][0] = atof(values[3]);
  angmom[m][1] = atof(values[4]);
  angmom[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecEM2::data_vel_hybrid(int m, char **values)
{
  angmom[m][0] = atof(values[0]);
  angmom[m][1] = atof(values[1]);
  angmom[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecEM2::pack_data(double **buf)
{
  double *quat;

  int m;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    m = 0;
    buf[i][m++] = ubuf(tag[i]).d;
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

int AtomVecEM2::pack_data_hybrid(int i, double *buf)
{
  double *quat = bonus[i].quat;

  int m = 0;
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

void AtomVecEM2::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %d %-1.16e %-1.16e %-1.16e %-1.16e "
            "%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4],buf[i][5],
            buf[i][6],buf[i][7],
            buf[i][8],buf[i][9],buf[i][10],
            (int) ubuf(buf[i][11]).i,(int) ubuf(buf[i][12]).i,
            (int) ubuf(buf[i][13]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecEM2::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e",
          buf[0],buf[1],buf[2],buf[3],buf[4],buf[5]);

  return 9;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecEM2::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = angmom[i][0];
    buf[i][5] = angmom[i][1];
    buf[i][6] = angmom[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecEM2::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = angmom[i][0];
  buf[1] = angmom[i][1];
  buf[2] = angmom[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecEM2::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT 
            " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecEM2::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ---------------------------------------------------------------------- */

void AtomVecEM2::force_clear(int n, size_t nbytes)
{
//  memset(&dphi[n][0],0,nbytes);
//  memset(&dphi[n][1],0,nbytes);
  memset(&dphi[n][0],0,2*nbytes);
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecEM2::property_atom(char *name)
{
  if (strcmp(name,"phim") == 0) return 0;
  if (strcmp(name,"phib") == 0) return 1;
  if (strcmp(name,"dphim") == 0) return 2;
  if (strcmp(name,"dphib") == 0) return 3;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecEM2::pack_property_atom(int index, double *buf,
                                     int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = phi[i][0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = phi[i][1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = dphi[i][0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = dphi[i][1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecEM2::memory_usage()
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

  bytes += nmax*sizeof(Bonus);

  return bytes;
}
