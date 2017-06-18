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

#include "compute_em2_energy.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "comm.h"
#include "modify.h"
#include "pair.h"
#include "force.h"
#include "string.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEM2Energy::ComputeEM2Energy(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute com command");

  int itmp;
  if (force->pair == NULL)
    error->all(FLERR,"compute em2_energy command is incompatible with Pair style");
  n_energy = *((int*) force->pair->extract("nenergy",itmp));
  energy = (double *) force->pair->extract("energy",itmp);
  if (!energy)
    error->all(FLERR,"compute em2_energy command is incompatible with Pair style");

  vector_flag = 1;
  size_vector = n_energy;
  extvector = 0;

  vector = new double[n_energy];
}

/* ---------------------------------------------------------------------- */

ComputeEM2Energy::~ComputeEM2Energy()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeEM2Energy::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_energy") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_energy");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2Energy::compute_vector()
{
  invoked_vector = update->ntimestep;

  for (int i=0;i<n_energy;++i) {
    vector[i] = energy[i];
  }
}
