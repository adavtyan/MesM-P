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

#ifdef COMPUTE_CLASS

ComputeStyle(em2_energy,ComputeEM2Energy)

#else

#ifndef LMP_COMPUTE_EM2_ENERGY_H
#define LMP_COMPUTE_EM2_ENERGY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEM2Energy : public Compute {
 public:
  ComputeEM2Energy(class LAMMPS *, int, char **);
  ~ComputeEM2Energy();
  void init();
  void compute_vector();

 private:
  int n_energy;
  double *energy;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
