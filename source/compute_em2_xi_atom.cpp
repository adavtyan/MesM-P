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
#include "compute_em2_xi_atom.h"
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

ComputeEM2XiAtom::ComputeEM2XiAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute em2_xhi/atom command");

  peratom_flag = 1;
  size_peratom_cols = 6;

  nmax = 0;
  xi_vec = NULL;
  xi_m = xi_b = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2XiAtom::~ComputeEM2XiAtom()
{
  memory->destroy(xi_vec);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2XiAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_xi/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_xi/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2XiAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow xi_vec array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(xi_vec);
    nmax = atom->nmax;
    memory->create(xi_vec,nmax,6,"xi_vec/atom:xi");
    array_atom = xi_vec;
  }

  // Extract xi_m and xi_b arrays from the Pair style
  if (!xi_m || !xi_b) {
   int itmp;
   if (force->pair == NULL)
     error->all(FLERR,"compute em2_phi/atom command is incompatible with Pair style");
   xi_m = (double **) force->pair->extract("xi_m",itmp);
   xi_b = (double **) force->pair->extract("xi_b",itmp);
   if (!xi_m || !xi_b)
     error->all(FLERR,"compute em2_phi/atom command is incompatible with Pair style");
  }


  // Pack xi array

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xi_vec[i][0] = xi_m[i][0];
      xi_vec[i][1] = xi_m[i][1];
      xi_vec[i][2] = xi_m[i][2];
      xi_vec[i][3] = xi_b[i][0];
      xi_vec[i][4] = xi_b[i][1];
      xi_vec[i][5] = xi_b[i][2];
    } else {
      xi_vec[i][0] = 0.0;
      xi_vec[i][1] = 0.0;
      xi_vec[i][2] = 0.0;
      xi_vec[i][3] = 0.0;
      xi_vec[i][4] = 0.0;
      xi_vec[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEM2XiAtom::memory_usage()
{
  double bytes = nmax * 6 * sizeof(double);
  return bytes;
}
