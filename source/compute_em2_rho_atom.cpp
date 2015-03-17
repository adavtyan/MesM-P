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

#include "string.h"
#include "compute_em2_rho_atom.h"
#include "atom.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEM2RhoAtom::ComputeEM2RhoAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute em2_rho/atom command");

  peratom_flag = 1;
  size_peratom_cols = 2;

  nmax = 0;
  rho_vec = NULL;
  rho_m = rho_b = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2RhoAtom::~ComputeEM2RhoAtom()
{
  memory->destroy(rho_vec);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2RhoAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_rho/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_rho/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2RhoAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow rho_vec array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(rho_vec);
    nmax = atom->nmax;
    memory->create(rho_vec,nmax,2,"rho_vec/atom:rho");
    array_atom = rho_vec;
  }

  // Extract rho_m and rho_b arrays from the Pair style
  if (!rho_m || !rho_b) {
    int itmp;
    if (force->pair == NULL)
      error->all(FLERR,"compute em2_phi/atom command is incompatible with Pair style");
    rho_m = (double *) force->pair->extract("rho_m",itmp);
    rho_b = (double *) force->pair->extract("rho_b",itmp);
    if (!rho_m || !rho_b)
      error->all(FLERR,"compute em2_phi/atom command is incompatible with Pair style");
  }

  // Pack rho array

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      rho_vec[i][0] = rho_m[i];
      rho_vec[i][1] = rho_b[i];
    } else {
      rho_vec[i][0] = 0.0;
      rho_vec[i][1] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEM2RhoAtom::memory_usage()
{
  double bytes = nmax * 2 * sizeof(double);
  return bytes;
}
