/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef FIX_CLASS

FixStyle(nvt/em2/rot,FixNVTEM2ROT)

#else

#ifndef LMP_FIX_NVT_EM2_ROT_H
#define LMP_FIX_NVT_EM2_ROT_H

#include "fix_nh_em2_rot.h"

namespace LAMMPS_NS {

class FixNVTEM2ROT : public FixNHEM2ROT {
 public:
  FixNVTEM2ROT(class LAMMPS *, int, char **);
  ~FixNVTEM2ROT();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control must be used with fix nvt/em2

Self-explanatory.

E: Pressure control can not be used with fix nvt/em2

Self-explanatory.

*/
