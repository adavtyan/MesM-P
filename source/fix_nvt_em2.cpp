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
#include "fix_nvt_em2.h"
#include "group.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVTEM2::FixNVTEM2(LAMMPS *lmp, int narg, char **arg) :
  FixNHEM2(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt/em2");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can not be used with fix nvt/em2");

  // create a new compute temp style
  // id = fix-ID + temp

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[4];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp/em2";
  newarg[3] = arg[narg-1];

  modify->add_compute(4,newarg);
  delete [] newarg;
  tcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

FixNVTEM2::~FixNVTEM2()
{
}
