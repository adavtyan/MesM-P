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
#include "compute_em2_nc_atom.h"
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

ComputeEM2NcAtom::ComputeEM2NcAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute em2_nc/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  nc_vec = NULL;
  nc = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2NcAtom::~ComputeEM2NcAtom()
{
  memory->destroy(nc_vec);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2NcAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_nc/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_nc/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2NcAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow nc_vec array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(nc_vec);
    nmax = atom->nmax;
    memory->create(nc_vec,nmax,"nc_vec/atom:nc");
    vector_atom = nc_vec;
  }

  // Extract n_c array from the Pair style
  if (!nc) {
    int itmp;
    if (force->pair == NULL)
      error->all(FLERR,"compute em2_nc/atom command is incompatible with Pair style");
    nc = (double *) force->pair->extract("n_c",itmp);
    if (!nc)
      error->all(FLERR,"compute em2_nc/atom command is incompatible with Pair style");
  }

  // Pack nc array

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      nc_vec[i] = nc[i];
    } else {
      nc_vec[i] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEM2NcAtom::memory_usage()
{
  double bytes = nmax * 2 * sizeof(double);
  return bytes;
}
