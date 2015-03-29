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

FixStyle(em2/spam,FixEM2SPAM)

#else

#ifndef LMP_FIX_EM2_SPAM_H
#define LMP_FIX_EM2_SPAM_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixEM2SPAM : public FixNVE {
 public:
  FixEM2SPAM(class LAMMPS *, int, char **);
  ~FixEM2SPAM();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void reset_dt();

 private:
  double dtv2;
  class AtomVecEM2 *avec;
  int allocated;
  char *conf_file;

  // per-type arrays
  double **inertia; // Not used here
  int *quat_flag; // Not used here
  int *mem_flag;
  int *prot_flag;

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
#endif

/* ERROR/WARNING messages:

E: Compute nve/em2 requires atom style em2

Self-explanatory.

*/
