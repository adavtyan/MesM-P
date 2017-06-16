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

FixStyle(bd,FixBD)

#else

#ifndef LMP_FIX_BD_H
#define LMP_FIX_BD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBD : public Fix {
 public:
  FixBD(class LAMMPS *, int, char **);
  virtual ~FixBD() {}
  int setmask();
  virtual void init();
  virtual void post_force(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void post_force_respa(int, int, int);
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void reset_dt();

 protected:
  double temp_bd, gamma;
  int seed;
  class RanPark *random;
  double dt;
  double *step_respa;
  int ilevel_respa;
  int mass_require;
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
