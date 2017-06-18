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

#ifndef LMP_FIX_NH_EM2_ROT_H
#define LMP_FIX_NH_EM2_ROT_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHEM2ROT : public FixNH {
 public:
  FixNHEM2ROT(class LAMMPS *, int, char **);
  virtual ~FixNHEM2ROT();
  void init();

 protected:
  class AtomVecEM2 *avec;
  int allocated;
  char *conf_file;

  void nve_v();
  void nve_x();
  void nh_v_temp();

  // per-type arrays
  double **inertia;
  int *quat_flag;
  int *mem_flag; // Not used here
  int *prot_flag; // Not used here

  void read_conf_file(char *filename);
  void allocate();
  char *ltrim(char *s);
  char *rtrim(char *s);
  char *trim(char *s);
  bool isEmptyString(char *str);
  inline void print_log(char *line);
};

}

#endif

/* ERROR/WARNING messages:

E: Compute nvt/nph/npt em2 requires atom style em2

Self-explanatory.

*/
