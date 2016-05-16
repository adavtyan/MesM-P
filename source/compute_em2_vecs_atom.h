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

#ifdef COMPUTE_CLASS

ComputeStyle(em2_vecs/atom,ComputeEM2VecsAtom)

#else

#ifndef LMP_COMPUTE_EM2_QUAT_VECS_H
#define LMP_COMPUTE_EM2_QUAT_VECS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEM2VecsAtom : public Compute {
 public:
  ComputeEM2VecsAtom(class LAMMPS *, int, char **);
  ~ComputeEM2VecsAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double **vecs;
  double scale_normal, scale_inplane;
  class AtomVecEM2 *avec;
};

}

#endif
#endif
