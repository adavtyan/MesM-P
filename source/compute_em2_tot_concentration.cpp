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

#include "mpi.h"
#include "compute_em2_tot_concentration.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "atom_vec_em2.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTotConc::ComputeTotConc(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute tot_conc command");

  vector_flag = 1;
  size_vector = 2;
  extscalar = 1;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeTotConc::~ComputeTotConc()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTotConc::init()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) avec = (AtomVecEM2 *) atom->style_match("em2_angle");
  if (!avec) error->all(FLERR,"compute em2_dphi/atom command requires atom style em2");
}

/* ---------------------------------------------------------------------- */

void ComputeTotConc::compute_vector()
{
  invoked_vector = update->ntimestep;

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **phi = avec->phi;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double sum[2] = {0.0, 0.0};

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      sum[0] -= 0.5*(phi[i][0] - 1.0);
      sum[1] -= 0.5*(phi[i][1] - 1.0);
    }
  }

  MPI_Allreduce(sum,vector,2,MPI_DOUBLE,MPI_SUM,world);
//  MPI_Allreduce(&sum[0],&vector[0],1,MPI_DOUBLE,MPI_SUM,world);
//  MPI_Allreduce(&sum[1],&vector[1],1,MPI_DOUBLE,MPI_SUM,world);
}
