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

ComputeStyle(em2_rho/atom,ComputeEM2RhoAtom)

#else

#ifndef LMP_COMPUTE_EM2_RHO_ATOM_H
#define LMP_COMPUTE_EM2_RHO_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEM2RhoAtom : public Compute {
 public:
  ComputeEM2RhoAtom(class LAMMPS *, int, char **);
  ~ComputeEM2RhoAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double **rho_vec;
  double *rho_m, *rho_b;
};

}

#endif
#endif
