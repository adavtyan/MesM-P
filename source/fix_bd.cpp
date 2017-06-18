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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "random_park.h"
#include "fix_bd.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBD::FixBD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix bd command");

  temp_bd = force->numeric(FLERR,arg[3]);
  gamma = force->numeric(FLERR,arg[4]);
  seed = force->inumeric(FLERR,arg[5]);

  random = new RanPark(lmp,seed + comm->me);

  dynamic_group_allow = 1;
  time_integrate = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

int FixBD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBD::init()
{
  dt = update->dt;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixBD::post_force(int vflag)
{
  // Sqrt(2*kT*gamma/dt) that are converted to force units (mass is included in gamma)
  // In the square root the kT is converted to the unites of mass*velocity^2 for consistancy
  double factor = sqrt(2.0*force->boltz*temp_bd*gamma/dt/force->mvv2e)/force->ftm2v;

  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal>0) random->reset(seed,x[0]);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += factor*random->gaussian();
      f[i][1] += factor*random->gaussian();
      f[i][2] += factor*random->gaussian();
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBD::initial_integrate(int vflag)
{
  // Nothing to be done here to BD
}

/* ---------------------------------------------------------------------- */

void FixBD::final_integrate()
{
  // dt/gamma that are converted from force*time/mass to velocity units
  double factor = force->ftm2v*dt/gamma;

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += f[i][0]*factor;
      x[i][1] += f[i][1]*factor;
      x[i][2] += f[i][2]*factor;
    }
  } 
}

/* ---------------------------------------------------------------------- */

void FixBD::post_force_respa(int vflag, int ilevel, int iloop)
{
  dt = step_respa[ilevel];

  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBD::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dt = step_respa[ilevel];

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixBD::final_integrate_respa(int ilevel, int iloop)
{
  dt = step_respa[ilevel];
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixBD::reset_dt()
{
  dt = update->dt;
}
