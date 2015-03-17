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
#include "compute_em2_phi_atom.h"
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

ComputeEM2PhiAtom::ComputeEM2PhiAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute em2_phi/atom command");

  peratom_flag = 1;
  size_peratom_cols = 2;

  nmax = 0;
  phi_vec = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2PhiAtom::~ComputeEM2PhiAtom()
{
  memory->destroy(phi_vec);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2PhiAtom::init()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) error->all(FLERR,"compute em2_phi/atom command requires atom style em2");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_phi/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_phi/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2PhiAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow phi_vec array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(phi_vec);
    nmax = atom->nmax;
    memory->create(phi_vec,nmax,2,"phi_vec/atom:phi");
    array_atom = phi_vec;
  }

  // Pack phi array

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **phi_half = avec->phi_half;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      phi_vec[i][0] = phi_half[i][0];
      phi_vec[i][1] = phi_half[i][1];
    } else {
      phi_vec[i][0] = 0.0;
      phi_vec[i][1] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEM2PhiAtom::memory_usage()
{
  double bytes = nmax * 2 * sizeof(double);
  return bytes;
}
