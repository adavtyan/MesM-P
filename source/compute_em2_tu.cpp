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

#include "compute_em2_tu.h"
#include "atom_vec_em2.h"
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

ComputeEM2TU::ComputeEM2TU(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute com command");

  scalar_flag = 1;
  extscalar = 0;

  nmax = 0;
  tu = NULL;
  nc = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2TU::~ComputeEM2TU()
{
  memory->destroy(tu);
  memory->destroy(nc);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2TU::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_tu") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_tu");

  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) error->all(FLERR,"Pair em2 requires atom style em2");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeEM2TU::compute_scalar()
{
/*  invoked_vector = update->ntimestep;

  int itmp;
  mem_stat = *((double*) force->pair->extract("mem_stat",itmp));
  prot_stat = *((double*) force->pair->extract("prot_stat",itmp));

  vector[0] = mem_stat;
  vector[1] = prot_stat;*/

  int i;

  if (atom->nmax > nmax) {
    memory->destroy(tu);
    memory->destroy(nc);
    nmax = atom->nmax;
    memory->create(tu,nmax,"compute:tu");
    memory->create(nc,nmax,"compute:nc");
  }

  for (i=0; i < nlocal; i++) {
    //
  }
}
