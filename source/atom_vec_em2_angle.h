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

#ifdef ATOM_CLASS

AtomStyle(em2_angle,AtomVecEM2Angle)

#else

#ifndef LMP_ATOM_VEC_EM2_ANGLE_H
#define LMP_ATOM_VEC_EM2_ANGLE_H

#include "atom_vec_em2.h"

namespace LAMMPS_NS {

class AtomVecEM2Angle : public AtomVecEM2 {
 public:

  AtomVecEM2Angle(class LAMMPS *);
  ~AtomVecEM2Angle() {}
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  int data_atom_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  bigint memory_usage();

 private:
  tagint *molecule;
  int **nspecial;
  tagint **special;
  int *num_bond;
  int **bond_type;
  tagint **bond_atom;
  int *num_angle;
  int **angle_type;
  tagint **angle_atom1,**angle_atom2,**angle_atom3;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

E: Assigning ellipsoid parameters to non-ellipsoid atom

Self-explanatory.

E: Invalid shape in Ellipsoids section of data file

Self-explanatory.

*/
