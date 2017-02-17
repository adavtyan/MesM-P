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
#include "compute_em2_dphi_atom.h"
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

ComputeEM2dPhiAtom::ComputeEM2dPhiAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute em2_dphi/atom command");

  peratom_flag = 1;
  size_peratom_cols = 2;

  nmax = 0;
  dphi_vec = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2dPhiAtom::~ComputeEM2dPhiAtom()
{
  memory->destroy(dphi_vec);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2dPhiAtom::init()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) avec = (AtomVecEM2 *) atom->style_match("em2_angle");
  if (!avec) error->all(FLERR,"compute em2_dphi/atom command requires atom style em2");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_dphi/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_dphi/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2dPhiAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow dphi_vec array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(dphi_vec);
    nmax = atom->nmax;
    memory->create(dphi_vec,nmax,2,"dphi_vec/atom:dphi");
    array_atom = dphi_vec;
  }

  // Pack dphi array

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **dphi = avec->dphi;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dphi_vec[i][0] = dphi[i][0];
      dphi_vec[i][1] = dphi[i][1];
    } else {
      dphi_vec[i][0] = 0.0;
      dphi_vec[i][1] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEM2dPhiAtom::memory_usage()
{
  double bytes = nmax * 2 * sizeof(double);
  return bytes;
}
