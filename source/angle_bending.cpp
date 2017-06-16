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

#include <math.h>
#include <stdlib.h>
#include "angle_bending.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleBending::AngleBending(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleBending::~AngleBending()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(kl);
    memory->destroy(lambda0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleBending::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double t1[3], t2[3], tdot;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,rdot,a1,a2;
  double r1r2_inv, rsq1_inv, rsq2_inv;
  double kl_r1r2_inv;

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond and normal vector t1

    delx1 = x[i2][0] - x[i1][0];
    dely1 = x[i2][1] - x[i1][1];
    delz1 = x[i2][2] - x[i1][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

/*    t1[0] = delx1/r1;
    t1[1] = dely1/r1;
    t1[2] = delz1/r1;*/

    // 2nd bond and normal vector t2

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

/*    t2[0] = delx2/r2;
    t2[1] = dely2/r2;
    t2[2] = delz2/r2;*/

    r1r2_inv = 1.0/(r1*r2);
    rsq1_inv = 1.0/rsq1;
    rsq2_inv = 1.0/rsq2;

    // dot product of r1 and r2 and t1 and t2
    rdot = delx1*delx2 + dely1*dely2 + delz1*delz2;
    tdot = rdot*r1r2_inv;
//    tdot = t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2];

    // force & energy    

    if (eflag) eangle = -kl[type]*tdot;

    a1 = rdot*rsq1_inv;
    a2 = rdot*rsq2_inv;

    kl_r1r2_inv = kl[type]*r1r2_inv;

    f1[0] = -kl_r1r2_inv*(delx2 - a1*delx1);
    f1[1] = -kl_r1r2_inv*(dely2 - a1*dely1);
    f1[2] = -kl_r1r2_inv*(delz2 - a1*delz1);

    f3[0] = -kl_r1r2_inv*(-delx1 + a2*delx2);
    f3[1] = -kl_r1r2_inv*(-dely1 + a2*dely2);
    f3[2] = -kl_r1r2_inv*(-delz1 + a2*delz2);

/*    f1[0] = kl[type]*(delx2 - a1*delx1)*r1r2_inv;
    f1[1] = kl[type]*(dely2 - a1*dely1)*r1r2_inv;
    f1[2] = kl[type]*(delz2 - a1*delz1)*r1r2_inv;

    f3[0] = kl[type]*(-delx1 + a2*delx2)*r1r2_inv;
    f3[1] = kl[type]*(-dely1 + a2*dely2)*r1r2_inv;
    f3[2] = kl[type]*(-delz1 + a2*delz2)*r1r2_inv;*/

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleBending::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(kl,n+1,"angle:kl");
  memory->create(lambda0,n+1,"angle:lambda0");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleBending::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double lambda0_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    lambda0[i] = lambda0_one;
    kl[i] = k_one/lambda0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleBending::equilibrium_angle(int i)
{
  return lambda0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleBending::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&lambda0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kl[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleBending::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nangletypes,fp);
    fread(&lambda0[1],sizeof(double),atom->nangletypes,fp);
    fread(&kl[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&lambda0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kl[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleBending::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],lambda0[i]);
}

/* ---------------------------------------------------------------------- */

double AngleBending::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i2][0] - x[i1][0];
  double dely1 = x[i2][1] - x[i1][1];
  double delz1 = x[i2][2] - x[i1][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double tdot = delx1*delx2 + dely1*dely2 + delz1*delz2;
  tdot /= r1*r2;

  return kl[type]*tdot;
}
