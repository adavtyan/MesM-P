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
#include "compute_em2_quat_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_em2.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEM2QuatAtom::ComputeEM2QuatAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute em2_quat/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;

  nmax = 0;
  quat_vec = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2QuatAtom::~ComputeEM2QuatAtom()
{
  memory->destroy(quat_vec);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2QuatAtom::init()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) error->all(FLERR,"compute em2_quat/atom command requires atom style em2");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_quat/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_quat/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2QuatAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow quat_vec array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(quat_vec);
    nmax = atom->nmax;
    memory->create(quat_vec,nmax,4,"quat_vec/atom:quat");
    array_atom = quat_vec;
  }

  // Pack quat array

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double *iquat;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      quat_vec[i][0] = bonus[i].quat[0];
      quat_vec[i][1] = bonus[i].quat[1];
      quat_vec[i][2] = bonus[i].quat[2];
      quat_vec[i][3] = bonus[i].quat[3];
    } else {
      quat_vec[i][0] = 0.0;
      quat_vec[i][1] = 0.0;
      quat_vec[i][2] = 0.0;
      quat_vec[i][3] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEM2QuatAtom::memory_usage()
{
  double bytes = nmax * 2 * sizeof(double);
  return bytes;
}
