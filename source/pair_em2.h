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

#ifdef PAIR_CLASS

PairStyle(em2,PairEM2)

#else

#ifndef LMP_PAIR_EM2_H
#define LMP_PAIR_EM2_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEM2 : public Pair {
 public:
  PairEM2(LAMMPS *lmp);
  virtual ~PairEM2();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void *extract(const char *, int &);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
//  void write_data(FILE *);
//  void write_data_all(FILE *);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 protected:
  enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

  double cut_global;
  double **cut;
  double norm0[3], nt0[3];

  // 2-16 LJ-like potential parameters
  double **lj216_pot_flag, **lj216_epsilon, **lj216_sigma, **lj216_k0;
  double **lj1,**lj2,**lj3,**lj4;

  // Lucy potential parameters
  int **lucy_pot_flag;
  double **lucy_epsilon, **lucy_sigma, **lucy_sigmasq;

  // Bending potential parameters
  int **bend_pot_flag;
  double **bend_epsilon, **bend_sigmasq, **bend_k0, **bend_gamma_epsilon;

  // Oligomerization potential parameters
  int **olig_pot_flag, **aolig;
  double **olig_epsilon;

  // Intrinsic Curvature Coupling potential parameters
  int *ic_pot_flag;
  double *ic_lambda_m, *ic_lambda_k, *ic_gamma_epsilon;

  // Composition Coupling potential parameters
  int *cc_pot_flag;
  double *cc_epsilon, *cc_zeta0;

  // Membrain Composition potential parameters
  int *mem_comp_pot_flag;
  double *mem_comp_xi_epsilon, *mem_comp_epsilon;
//  double mem_comp_rcut_inv, mem_comp_rcutsq;

  // Protein Composition potential parameters
  int *prot_comp_pot_flag;
  double *prot_comp_xi_epsilon, *prot_comp_epsilon;
//  double prot_comp_rcut_inv, prot_comp_rcutsq;

  // Global parameters for Composition Potentials
  int comp_flag;
  double comp_rcut, comp_rcut_inv, comp_rcutsq;

  // Composition stat_variables
  int *mem_stat_flag, *prot_stat_flag;
  double mem_stat_epsilon, prot_stat_epsilon;
  double mem_stat, prot_stat; 

  // SPAM variables
  int spam_flag; 
  int flow_term_flag;
  double spam_gamma;

  // Variables for density and composition gradient calculations
  int nmax;
  double *rho_m, *rho_b, **xi_m, **xi_b; // Per-atom arrays
  int comm_ind; // Index of communication procedure to invoke

  double **offset;
  class AtomVecEM2 *avec;

  char *parfile;
  bool read_par_flag;

  void allocate();
  void read_parameters();
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

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair gayberne requires atom style ellipsoid

Self-explanatory.

E: Pair gayberne requires atoms with same type have same shape

Self-explanatory.

E: Pair gayberne epsilon a,b,c coeffs are not all set

Each atom type involved in pair_style gayberne must
have these 3 coefficients set at least once.

E: Bad matrix inversion in mldivide3

This error should not occur unless the matrix is badly formed.

*/
