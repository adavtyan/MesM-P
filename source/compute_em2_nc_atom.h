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

ComputeStyle(em2_nc/atom,ComputeEM2NcAtom)

#else

#ifndef LMP_COMPUTE_EM2_NC_ATOM_H
#define LMP_COMPUTE_EM2_NC_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEM2NcAtom : public Compute {
 public:
  ComputeEM2NcAtom(class LAMMPS *, int, char **);
  ~ComputeEM2NcAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double *nc_vec;
  double *nc;
};

}

#endif
#endif
