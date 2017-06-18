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

#include "compute_em2_stat.h"
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

ComputeEM2Stat::ComputeEM2Stat(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute com command");

  if (force->pair == NULL)
    error->all(FLERR,"compute em2_stat command is incompatible with Pair style");

  if (strcmp(force->pair_style,"em2") != 0)
    error->all(FLERR,"compute em2_stat command is incompatible with Pair style");

  vector_flag = 1;
  size_vector = 2;
  extvector = 0;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeEM2Stat::~ComputeEM2Stat()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeEM2Stat::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_stat") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_stat");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2Stat::compute_vector()
{
  invoked_vector = update->ntimestep;

  int itmp;
  mem_stat = *((double*) force->pair->extract("mem_stat",itmp));
  prot_stat = *((double*) force->pair->extract("prot_stat",itmp));

  vector[0] = mem_stat;
  vector[1] = prot_stat;
}
