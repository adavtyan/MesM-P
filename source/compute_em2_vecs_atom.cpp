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
   Description: compute em2_vecs return starting and ending points of 
                two vectors that point in the direction of normals and 
                in-plain vectors of the membrane quasi-particles.
                Each vector can optionally be scaled by compute command
                argumnets. If one scaling argumnet is given, both vectors 
                are scaled equally, while if two scaling argumnets are 
                given, the normal and in-plane vectors are scaled using them,
                accordingely.
   Examples: compute    id mem_group em2_vecs/atom
             compute    id mem_group em2_vecs/atom 2.0 
             compute    id mem_group em2_vecs/atom 2.0 4.0 
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_em2_vecs_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_em2.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEM2VecsAtom::ComputeEM2VecsAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3 || narg>5) error->all(FLERR,"Illegal compute em2_vecs/atom command");

  peratom_flag = 1;
  size_peratom_cols = 12;

  scale_normal = scale_inplane = 1.0;

  if (narg==4) {
    scale_normal = scale_inplane = force->numeric(FLERR,arg[3]);
  } else if (narg==5) {
    scale_normal = force->numeric(FLERR,arg[3]);
    scale_inplane = force->numeric(FLERR,arg[4]);
  }

  nmax = 0;
  vecs = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEM2VecsAtom::~ComputeEM2VecsAtom()
{
  memory->destroy(vecs);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2VecsAtom::init()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) avec = (AtomVecEM2 *) atom->style_match("em2_angle");
  if (!avec) error->all(FLERR,"compute em2_vecs/atom command requires atom style em2");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_vecs/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_vecs/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEM2VecsAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow vecs array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(vecs);
    nmax = atom->nmax;
    memory->create(vecs,nmax,12,"em2_vecs/atom:vecs");
    array_atom = vecs;
  }

  // Compute and pack vecs array

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double *iquat;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double ai[3][3],normi[3],nti[3];
  double norm0[3] = {0.0, 0.0, 1.0};
  double nt0[3] = {0.0, 1.0, 0.0};

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      iquat = bonus[i].quat;
      MathExtra::quat_to_mat_trans(iquat, ai);
      MathExtra::vecmat(norm0,ai,normi);
      MathExtra::vecmat(nt0,ai,nti);

      vecs[i][0] = x[i][0];
      vecs[i][1] = x[i][1];
      vecs[i][2] = x[i][2];

      vecs[i][3] = vecs[i][0] + scale_normal*normi[0];
      vecs[i][4] = vecs[i][1] + scale_normal*normi[1];
      vecs[i][5] = vecs[i][2] + scale_normal*normi[2];

      vecs[i][6] = x[i][0];
      vecs[i][7] = x[i][1];
      vecs[i][8] = x[i][2];

      vecs[i][9] = vecs[i][6] + scale_inplane*nti[0];
      vecs[i][10] = vecs[i][7] + scale_inplane*nti[1];
      vecs[i][11] = vecs[i][8] + scale_inplane*nti[2];
    } else {
      vecs[i][0] = vecs[i][1] = vecs[i][2] = 0.0;
      vecs[i][3] = vecs[i][4] = vecs[i][5] = 0.0;
      vecs[i][6] = vecs[i][7] = vecs[i][8] = 0.0;
      vecs[i][9] = vecs[i][10] = vecs[i][11] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEM2VecsAtom::memory_usage()
{
  double bytes = nmax * 2 * sizeof(double);
  return bytes;
}
