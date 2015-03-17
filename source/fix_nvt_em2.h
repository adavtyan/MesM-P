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

#ifdef FIX_CLASS

FixStyle(nvt/em2,FixNVTEM2)

#else

#ifndef LMP_FIX_NVT_EM2_H
#define LMP_FIX_NVT_EM2_H

#include "fix_nh_em2.h"

namespace LAMMPS_NS {

class FixNVTEM2 : public FixNHEM2 {
 public:
  FixNVTEM2(class LAMMPS *, int, char **);
  ~FixNVTEM2();
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
