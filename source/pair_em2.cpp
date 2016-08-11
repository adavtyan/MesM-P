/* ----------------------------------------------------------------------
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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_em2.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_em2.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"

#include "update.h"

using namespace LAMMPS_NS;

/*static const char cite_pair_em2[] =
  "pair em2 command:\n\n"
  "@Article{AuthorYY,\n"
  " author =  {A. A. Author1, A. A. Author2, ...},\n"
  " title =   {Tittle},\n"
  " journal = {Journal e.g. J.~Chem.~Phys.},\n"
  " year =    YYYY,\n"
  " volume =  vol,\n"
  " pages =   {pages}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairEM2::PairEM2(LAMMPS *lmp) : Pair(lmp)
{ 
//  if (lmp->citeme) lmp->citeme->add(cite_pair_em2);

  single_enable = 0;
  writedata = 0;

  nmax = 0;
  rho_m = NULL;
  rho_b = NULL;
  xi_m = NULL;
  xi_b = NULL;
  n_c = NULL;

  mem_comp_poly_exp = NULL;
  mem_comp_poly_coeff = NULL;

  // set comm size needed by this Pair

  comm_forward = 6;
  comm_reverse = 6;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairEM2::~PairEM2()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(offset);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(lj216_pot_flag);
    memory->destroy(lj216_epsilon);
    memory->destroy(lj216_sigma);
    memory->destroy(lj216_k0);
    memory->destroy(lj612_1);
    memory->destroy(lj612_2);
    memory->destroy(lj612_3);
    memory->destroy(lj612_4);
    memory->destroy(lj612_pot_flag);
    memory->destroy(lj612_epsilon);
    memory->destroy(lj612_sigma);
    memory->destroy(lj612_k0);
    memory->destroy(lucy_pot_flag);
    memory->destroy(lucy_epsilon);
    memory->destroy(lucy_sigma);
    memory->destroy(lucy_sigma_inv);
    memory->destroy(lucy_sigmasq_inv);
    memory->destroy(lucy_table_pot_flag);
    memory->destroy(lucy_table_epsilon);
    memory->destroy(lucy_table_sigma);
    memory->destroy(ex12_pot_flag);
    memory->destroy(ex12_epsilon);
    memory->destroy(ex12_sigma);
    memory->destroy(ex12_epsig);
    memory->destroy(gauss_pot_flag);
    memory->destroy(gauss_epsilon);
    memory->destroy(gauss_sigma);
    memory->destroy(gauss_sigmasq_inv);
    memory->destroy(gauss_r0);
    memory->destroy(bend_pot_flag);
    memory->destroy(bend_epsilon);
    memory->destroy(bend_k0);
    memory->destroy(bend_gamma_epsilon);
    memory->destroy(lipid_factor_flag);
    memory->destroy(lipid_lambda);
    memory->destroy(olig_pot_flag);
    memory->destroy(olig_epsilon);
    memory->destroy(aolig);
    memory->destroy(ic_pot_flag);
    memory->destroy(ic_lambda_m);
    memory->destroy(ic_lambda_k);
    memory->destroy(ic_gamma_epsilon);
    memory->destroy(cc_pot_flag);
    memory->destroy(cc_epsilon);
    memory->destroy(cc_zeta0);
    memory->destroy(mem_comp_pot_flag);
    memory->destroy(mem_comp_epsilon);
    memory->destroy(mem_comp_xi_epsilon);
    memory->destroy(mem_comp_npoly);
    memory->destroy(prot_comp_pot_flag);
    memory->destroy(prot_comp_epsilon);
    memory->destroy(prot_comp_xi_epsilon);
    memory->destroy(mem_stat_flag);
    memory->destroy(prot_stat_flag);
    memory->destroy(pcov_pot_flag);

    memory->destroy(rho_m);
    memory->destroy(rho_b);
    memory->destroy(xi_m);
    memory->destroy(xi_b);
    memory->destroy(n_c);

    int i, j;

    for (i = 1; i <= atom->ntypes; i++) {
      if (mem_comp_poly_exp[i]!=NULL) delete [] mem_comp_poly_exp[i];
      if (mem_comp_poly_coeff[i]!=NULL) delete [] mem_comp_poly_coeff[i];
    }
    if (mem_comp_poly_exp!=NULL) delete [] mem_comp_poly_exp;
    if (mem_comp_poly_coeff!=NULL) delete [] mem_comp_poly_coeff;

    for (int i = 1; i <= atom->ntypes; i++) {
      for (int j = i; j <= atom->ntypes; j++) {
        if (lucy_etable1[i][j]) free(lucy_etable1[i][j]);
        if (lucy_etable2[i][j]) free(lucy_etable2[i][j]);
        if (lucy_ftable1[i][j]) free(lucy_ftable1[i][j]);
        if (lucy_ftable2[i][j]) free(lucy_ftable2[i][j]);
      }
      free(lucy_etable1[i]);
      free(lucy_etable2[i]);
      free(lucy_ftable1[i]);
      free(lucy_ftable2[i]);
    }
    free(lucy_etable1);
    free(lucy_etable2);
    free(lucy_ftable1);
    free(lucy_ftable2);
  }
}

/* ---------------------------------------------------------------------- */

void PairEM2::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl,one_eng,r,rsq,rinv,r2inv,r4inv,r6inv,r12inv,r16inv,fpair,factor_lj,epsilon,sigma,eng,eng_lj;
  double fforce[3],gb_force[3],torque[3],ttor[3],rtor[3],r12[3],ru12[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp;

  // Quaternions and orientations
  double *iquat,*jquat;
  double ai[3][3],aj[3][3],normi[3],normj[3],nti[3],ntj[3];

  // SPAM
  double K_factor, K_force_factor, lipid_factor;
  double fij, fij_force_factor;
  bool b_K_force, b_fij_force;
  bool b_K_force_i, b_K_force_j;
  bool b_fij_force_i, b_fij_force_j;
  bool b_spam_force, b_spam_force2, m_spam_force;
  bool b_spam_force_i, b_spam_force_j;
  double iphi_m,jphi_m,iphi_b,jphi_b;
  double dpm,dpb,idphi[2],jdphi[2];
  double rho_one,rho_ave_inv,xi_one,dw;
  // For Stat calculations
  double sum_phi_m, sum_phi_b, sum_phi_m_all, sum_phi_b_all;
  double n_phi_m, n_phi_b, n_phi_m_all, n_phi_b_all;
  // For Protein Coverage potential
  double pcov_sum_phi_b, pcov_sum_phi_b_all;
  double pcov_n_phi_b, pcov_n_phi_b_all;
  double pcov_n_sol_phi_b, pcov_n_sol_phi_b_all;

  // Lucy potential
  double epsq2;

  // Lucy table potential
  int ir;
  double *etb1, *etb2, *ftb1, *ftb2;

  // Bending potential
  double uij, gamma,gamma_r, gamma_epsilon, thetai, thetaj, alphai, alphaj, alpha_sq;
  double epsilon_sigmasq_r2inv, thetai_gamma_r, thetaj_gamma_r, fi, fj;
  double fthi, fthj, fali, falj, tor_thi, tor_thj, tor_ali, tor_alj;

  // Composition coupling potential
  double zeta_m, zeta_b;

  // Composition potentials
  double phi2, phi4;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **phi = avec->phi;
  double **dphi = avec->dphi;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // Number of atoms involved in pair calculations
  int npall = newton_pair ? nall : nlocal;

  for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // calculate per-atom densities and gradients
  // if Membrain Composition or Protein Composition potentials are on

//  if (spam_flag && comp_flag) {

    // grow energy array if necessary

    if (atom->nmax > nmax) {
      memory->destroy(rho_m);
      memory->destroy(rho_b);
      memory->destroy(xi_m);
      memory->destroy(xi_b);
      memory->destroy(n_c);
      nmax = atom->nmax;
      memory->create(rho_m,nmax,"pair:rho_m");
      memory->create(rho_b,nmax,"pair:rho_b");
      memory->create(xi_m,nmax,3,"pair:xi_m");
      memory->create(xi_b,nmax,3,"pair:xi_b");
      memory->create(n_c,nmax,"pair:n_c");
    }

    // zero out per-atom density, gradient, and nuber of neighbours

    for (i = 0; i < npall; i++) {
      rho_m[i] = 0.0;
      rho_b[i] = 0.0;
      xi_m[i][0] = 0.0;
      xi_m[i][1] = 0.0;
      xi_m[i][2] = 0.0;
      xi_b[i][0] = 0.0;
      xi_b[i][1] = 0.0;
      xi_b[i][2] = 0.0;
      n_c[i] = 0.0;
    }

    // calculate densities rho_m and rho_b

    // First calculating self term: m*W(0) = m
    for (i = 0; i < nlocal; i++) {
      itype = type[i];
      if (mem_comp_pot_flag[itype]) rho_m[i] = mass[itype];
      if (prot_comp_pot_flag[itype]) rho_b[i] = mass[itype];
    }

    // loop over neighbors of my atoms

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        // r12 = center to center vector

        r12[0] = xtmp - x[j][0];
        r12[1] = ytmp - x[j][1];
        r12[2] = ztmp - x[j][2];
        rsq = MathExtra::dot3(r12,r12);
        jtype = type[j];

        // compute if less than cutoff

        // Nummber of neighbours calculation
        if (rsq < cutsq[itype][jtype] &&  ic_pot_flag[itype] && ic_pot_flag[jtype]) {
          n_c[i] += 1.0;
          if (newton_pair || j < nlocal)
            n_c[j] += 1.0;
        }

        // Density calculation
        if (rsq < comp_rcutsq) {
          r = sqrt(rsq);
          sigma = r*comp_rcut_inv;
          rho_one =(1.0-sigma)*(1.0-sigma)*(1.0-sigma)*(1.0+3.0*sigma);

          // Calculate density of particles carrying membrain composition
          if (mem_comp_pot_flag[itype] && mem_comp_pot_flag[jtype]) {
            rho_m[i] += mass[jtype]*rho_one;
            if (newton_pair || j < nlocal) {
              rho_m[j] += mass[itype]*rho_one;
            }
          }

          // Calculate density of particles carrying protein composition
          if (prot_comp_pot_flag[itype] && prot_comp_pot_flag[jtype]) {
            rho_b[i] += mass[jtype]*rho_one;
            if (newton_pair || j < nlocal) {
              rho_b[j] += mass[itype]*rho_one;
            }
          }
        }
      }
    }

    // communicate and sum densities

    comm_ind = 1;
    if (newton_pair) comm->reverse_comm_pair(this);
    comm->forward_comm_pair(this);

    // calculate gradients xi_m and xi_b
    // loop over neighbors of my atoms

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // SPAM variables, phi array
      // 0 - membrain, 1 - protein
      iphi_m = phi[i][0];
      iphi_b = phi[i][1];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        // r12 = center to center vector

        r12[0] = xtmp - x[j][0];
        r12[1] = ytmp - x[j][1];
        r12[2] = ztmp - x[j][2];
        rsq = MathExtra::dot3(r12,r12);
        jtype = type[j];

        // compute if less than cutoff

        // Composition gradient calulation
        if (rsq < comp_rcutsq) {

          // SPAM variables, phi array
          // 0 - membrain, 1 - protein
          jphi_m = phi[j][0];
          jphi_b = phi[j][1];

          r = sqrt(rsq);
          sigma = r*comp_rcut_inv; 
          dw = -12.0*(1.0 - sigma)*(1.0 - sigma)*comp_rcut_inv*comp_rcut_inv;

          // calculating gradient for membrain
          // Unlike the older EM2 code, rho_m (rather then rho_b) is used for xi_m calculations here

          if (mem_comp_pot_flag[itype] && mem_comp_pot_flag[jtype]) {
            rho_ave_inv = 1.0/(0.5*(rho_m[i]+rho_m[j]));
            xi_one = rho_ave_inv*(jphi_m - iphi_m)*dw;

            xi_m[i][0] += mass[jtype]*xi_one*r12[0];
            xi_m[i][1] += mass[jtype]*xi_one*r12[1];
            xi_m[i][2] += mass[jtype]*xi_one*r12[2];
            if (newton_pair || j < nlocal) {
              xi_m[j][0] += mass[itype]*xi_one*r12[0];
              xi_m[j][1] += mass[itype]*xi_one*r12[1];
              xi_m[j][2] += mass[itype]*xi_one*r12[2];
            }
          }

          // calculating gradient for protein

          if (prot_comp_pot_flag[itype] && prot_comp_pot_flag[jtype]) {
            rho_ave_inv = 1.0/(0.5*(rho_b[i]+rho_b[j]));
            xi_one = rho_ave_inv*(jphi_b - iphi_b)*dw;

            xi_b[i][0] += mass[jtype]*xi_one*r12[0];
            xi_b[i][1] += mass[jtype]*xi_one*r12[1];
            xi_b[i][2] += mass[jtype]*xi_one*r12[2];
            if (newton_pair || j < nlocal) {
              xi_b[j][0] += mass[itype]*xi_one*r12[0];
              xi_b[j][1] += mass[itype]*xi_one*r12[1];
              xi_b[j][2] += mass[itype]*xi_one*r12[2];
            }
          }
        }
      }
    }

    // communicate and sum gradients

    comm_ind = 2;
    if (newton_pair) comm->reverse_comm_pair(this);
    comm->forward_comm_pair(this);
//  }

  // Main loop
  // loop over neighbors of my atoms
  // idphi and jdphi are multipied by spam_gamma in the end

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    // Normal and in-plain vectors for particle i
//    if (bend_pot_flag[itype][jtype] || olig_pot_flag[itype][jtype]) {
      iquat = bonus[i].quat;
      MathExtra::quat_to_mat_trans(iquat, ai);
      MathExtra::vecmat(norm0,ai,normi);
      MathExtra::vecmat(nt0,ai,nti);
//    }

    // SPAM variables, phi array
    // 0 - membrain, 1 - protein
    if (spam_flag) {
      iphi_m = phi[i][0];
      iphi_b = phi[i][1];
    }

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // r12 = center to center vector

      r12[0] = xtmp - x[j][0];
      r12[1] = ytmp - x[j][1];
      r12[2] = ztmp - x[j][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {

        // Normal and in-plain vectors for particle j
//        if (bend_pot_flag[itype][jtype] || olig_pot_flag[itype][jtype]) {
          jquat = bonus[j].quat;
          MathExtra::quat_to_mat_trans(jquat, aj);
          MathExtra::vecmat(norm0,aj,normj);
          MathExtra::vecmat(nt0,aj,ntj);
//        }

        // SPAM variables, phi array
        // 0 - membrain, 1 - protein
        if (spam_flag) {
          jphi_m = phi[j][0];
          jphi_b = phi[j][1];
        }

        one_eng = 0.0;
        fforce[0] = fforce[1] = fforce[2] = 0.0;
        ttor[0] = ttor[1] = ttor[2] = 0.0;
        rtor[0] = rtor[1] = rtor[2] = 0.0;
        idphi[0] = idphi[1] = 0.0;
        jdphi[0] = jdphi[1] = 0.0;

        r = sqrt(rsq);

        // Common spam pre-factor used by most potentials, and common boolean functions
        // When Kpd_flag or Gpd_flag is 1, the full range (i.e. [-1,1]) of protein composition is used
        // When kpd_flag or Gpd_flag is 0, the potentials are varied only in [-1,0] range of protein composition
        if (Kpd_flag) {
          K_factor = 0.25*(MIN(MAX(iphi_b,-1.0),1.0) + MIN(MAX(jphi_b,-1.0),1.0) - 2.0);
          K_force_factor = 0.25;
          b_K_force_i = iphi_b<=1.0 && iphi_b>=-1.0;
          b_K_force_j = jphi_b<=1.0 && jphi_b>=-1.0;
        } else {
          K_factor = 0.5*(MIN(MAX(iphi_b,-1.0),0.0) + MIN(MAX(jphi_b,-1.0),0.0));
          K_force_factor = 0.5;
          b_K_force_i = iphi_b<=0.0 && iphi_b>=-1.0;
          b_K_force_j = jphi_b<=0.0 && jphi_b>=-1.0;
        }
        if (Gpd_flag) {
          fij = 0.125*(MAX(MIN(-iphi_b,1.0),-1.0) + MAX(MIN(-jphi_b,1.0),-1.0) + 2.0);
          fij_force_factor = 0.125;
          b_fij_force_i = iphi_b<=1.0 && iphi_b>=-1.0;
          b_fij_force_j = jphi_b<=1.0 && jphi_b>=-1.0;
        } else {
          fij = 0.25*(MAX(MIN(-iphi_b,1.0),0.0) + MAX(MIN(-jphi_b,1.0),0.0));
          fij_force_factor = 0.25;
          b_fij_force_i = iphi_b<=0.0 && iphi_b>=-1.0;
          b_fij_force_j = jphi_b<=0.0 && jphi_b>=-1.0;
        }
        b_K_force = b_K_force_i || (newton_pair && b_K_force_j);
        b_fij_force = b_fij_force_i || (newton_pair && b_fij_force_j);
        b_spam_force_i = b_K_force_i || b_fij_force_i;
        b_spam_force_j = b_K_force_j || b_fij_force_j;
        b_spam_force = b_K_force || b_fij_force;
        b_spam_force2 = (iphi_b<=1.0 && iphi_b>=-1.0) || (newton_pair && jphi_b<=1.0 && jphi_b>=-1.0);

        // 2-16 LJ-like term for membrane-membrane interaction
        // Forces, energy, and dphi_b are calculated
        if (lj216_pot_flag[itype][jtype]) {
          epsilon = 1.0;
          if (spam_flag) epsilon *= (1 - lj216_k0[itype][jtype]*K_factor);

          r2inv = 1.0/rsq;
          r4inv = r2inv*r2inv;
          r16inv = r4inv*r4inv*r4inv*r4inv;
          fpair = epsilon*(r16inv*lj1[itype][jtype]-r2inv*lj2[itype][jtype]);
          fpair *= r2inv;

          eng_lj = (r16inv*lj3[itype][jtype]-r2inv*lj4[itype][jtype]);
          eng = epsilon*eng_lj;

          energy[ET_LJ216] += eng;
          if (eflag) one_eng += eng;
          fforce[0] += r12[0]*fpair;
          fforce[1] += r12[1]*fpair;
          fforce[2] += r12[2]*fpair;

          // SPAM variable dynamics
          // Calculating -dphi/dt for phi_b only
          if (spam_flag && b_K_force) {
            // Contribution from the derivative of prefactor
            dpb = K_force_factor*lj216_k0[itype][jtype]*eng_lj;

            if (b_K_force_i) idphi[1] += dpb;
            if (b_K_force_j) jdphi[1] += dpb;
          }
        }

        // 6-12 LJ-like term for membrane-membrane interaction
        // Forces, energy, and dphi_b are calculated
        if (lj612_pot_flag[itype][jtype]) {
          epsilon = 1.0;
          if (spam_flag) epsilon *= (1 - lj612_k0[itype][jtype]*K_factor);

          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          r12inv = r6inv*r6inv;
          fpair = epsilon*(r12inv*lj612_1[itype][jtype]-r6inv*lj612_2[itype][jtype]);
          fpair *= r2inv;

          eng_lj = (r12inv*lj612_3[itype][jtype]-r6inv*lj612_4[itype][jtype]);
          eng = epsilon*eng_lj;

          energy[ET_LJ612] += eng;
          if (eflag) one_eng += eng;
          fforce[0] += r12[0]*fpair;
          fforce[1] += r12[1]*fpair;
          fforce[2] += r12[2]*fpair;

          // SPAM variable dynamics
          // Calculating -dphi/dt for phi_b only
          if (spam_flag && b_K_force) {
            // Contribution from the derivative of prefactor
            dpb = K_force_factor*lj612_k0[itype][jtype]*eng_lj;

            if (b_K_force_i) idphi[1] += dpb;
            if (b_K_force_j) jdphi[1] += dpb;
          }
        }

        // Repulsive term for solvent-solvent and membrane-solvet interaction
        // This term has form of Lucy kernel function
        // Only forces and energy are calculated
        if (lucy_pot_flag[itype][jtype] && r<lucy_sigma[itype][jtype]) {
          epsilon = lucy_epsilon[itype][jtype];

          sigma = r*lucy_sigma_inv[itype][jtype];
          epsq2 = epsilon*(1-sigma)*(1-sigma);
          fpair = 12*epsq2*lucy_sigmasq_inv[itype][jtype];

          eng = epsq2*(1-sigma)*(1+3*sigma);

          energy[ET_LUCY] += eng;
          if (eflag) one_eng += eng;
          fforce[0] += r12[0]*fpair;
          fforce[1] += r12[1]*fpair;
          fforce[2] += r12[2]*fpair; 
        }

        // Table version of Lucy potential
        // Repulsive term for solvent-solvent and membrane-solvet interaction
        // This term has form of Lucy kernel function
        // Only forces and energy are calculated
        if (lucy_table_pot_flag[itype][jtype] && r<lucy_table_sigma[itype][jtype]) {
          etb1 = lucy_etable1[itype][jtype];
          etb2 = lucy_etable2[itype][jtype];
          ftb1 = lucy_ftable1[itype][jtype];
          ftb2 = lucy_ftable2[itype][jtype];
          ir = int(r*lucy_tb_dr_inv);

          // Energy and force values are obtained from trangle interpolation
          
          eng = etb1[ir]*r + etb2[ir];
          fpair = ftb1[ir]*r + ftb2[ir];

          energy[ET_LUCY] += eng;
          if (eflag) one_eng += eng;
          fforce[0] += r12[0]*fpair;
          fforce[1] += r12[1]*fpair;
          fforce[2] += r12[2]*fpair; 
        }

        // Repulsive 1/r^12 term
        // Only forces and energy are calculated
        if (ex12_pot_flag[itype][jtype]) {
          epsilon = ex12_epsig[itype][jtype];

          r2inv = 1.0/rsq;
          r4inv = r2inv*r2inv;
          r12inv = r4inv*r4inv*r4inv;
          fpair = 12.0*epsilon*r12inv;
          fpair *= r2inv;

          eng = epsilon*r12inv;

          energy[ET_EX12] += eng;
          if (eflag) one_eng += eng;
          fforce[0] += r12[0]*fpair;
          fforce[1] += r12[1]*fpair;
          fforce[2] += r12[2]*fpair;
        }

        // Attractive gaussian potential
        // Only forces and energy are calculated
        if (gauss_pot_flag[itype][jtype]) {
          epsilon = gauss_epsilon[itype][jtype];

          eng = -epsilon*exp(-0.5*gauss_sigmasq_inv[itype][jtype]*(r-gauss_r0[itype][jtype])*(r-gauss_r0[itype][jtype]));
          fpair = eng*gauss_sigmasq_inv[itype][jtype]*(1.0 - gauss_r0[itype][jtype]/r);

          energy[ET_GAUSS] += eng;
          if (eflag) one_eng += eng;
          fforce[0] += r12[0]*fpair;
          fforce[1] += r12[1]*fpair;
          fforce[2] += r12[2]*fpair;
        }

        // Membrane elastic bending energy term
        // Forces, torques, energy, and dphi_b are calculated
        if (bend_pot_flag[itype][jtype]) {

          rinv = 1.0/r;
          r2inv = 1.0/rsq;

          // 4.0*epsilon*(1 - k0*( MAX(iphi_b,1.0) + MAX(jphi_b,1.0) ))*(sigma/r)^2
          epsilon_sigmasq_r2inv = bend_epsilon[itype][jtype]*r2inv;
          if (spam_flag) epsilon_sigmasq_r2inv *= (1 - bend_k0[itype][jtype]*K_factor);

          // Multiply by lipid rigidity bending factor
          // 0.5*(1 + lambda) + 0.5*(lampda - 1) * 0.5*(phi_m_i + phi_m_j)
          if (spam_flag && lipid_factor_flag[itype][jtype]) {
            lipid_factor = 0.5*(MIN(MAX(iphi_m,-1.0),1.0) + MIN(MAX(jphi_m,-1.0),1.0));
            lipid_factor = 0.5*(1.0 + lipid_lambda[itype][jtype] + lipid_factor * (lipid_lambda[itype][jtype] - 1.0));
            m_spam_force = (iphi_m<=1.0 && iphi_m>=-1.0) || (newton_pair && jphi_m<=1.0 && jphi_m>=-1.0);

            epsilon_sigmasq_r2inv *= lipid_factor;
          } else {
            lipid_factor = 1.0;
            m_spam_force = false;
          }

          ru12[0] = r12[0]*rinv;
          ru12[1] = r12[1]*rinv;
          ru12[2] = r12[2]*rinv;

          thetai = MathExtra::dot3(normi,ru12);
          thetaj = MathExtra::dot3(normj,ru12);

          alphai = MathExtra::dot3(nti,ru12);
          alphaj = MathExtra::dot3(ntj,ru12);

          alpha_sq = alphai*alphai + alphaj*alphaj;

          gamma_epsilon  = 0.5*bend_gamma_epsilon[itype][jtype];
          if (spam_flag) {
            //fi = MAX(MIN(-iphi_b,1.0),0.0);
            //fj = MAX(MIN(-jphi_b,1.0),0.0);
            //gamma_epsilon *= 0.25*(fi + fj);
            gamma_epsilon *= fij;
          }
          gamma = gamma_epsilon*alpha_sq;
          gamma_r = gamma*r;

          thetai_gamma_r = thetai - gamma_r;
          thetaj_gamma_r = thetaj + gamma_r;

          uij = thetai_gamma_r*thetai_gamma_r + thetaj_gamma_r*thetaj_gamma_r;

          // Spherical force
          fpair = (uij + thetai_gamma_r*thetai + thetaj_gamma_r*thetaj)*r2inv;
          fpair += (thetai_gamma_r - thetaj_gamma_r)*(gamma - 2.0*gamma_epsilon*alpha_sq)*rinv;
          fpair *= 2.0*epsilon_sigmasq_r2inv;

          fthi = -thetai_gamma_r*rinv;
          fthj = -thetaj_gamma_r*rinv;

          fali = 2.0*gamma_epsilon*(thetai_gamma_r - thetaj_gamma_r)*alphai;
          falj = 2.0*gamma_epsilon*(thetai_gamma_r - thetaj_gamma_r)*alphaj;

          // Aspherical forces
          gb_force[0] = fthi*normi[0] + fthj*normj[0] + fali*nti[0] + falj*ntj[0];
          gb_force[1] = fthi*normi[1] + fthj*normj[1] + fali*nti[1] + falj*ntj[1];
          gb_force[2] = fthi*normi[2] + fthj*normj[2] + fali*nti[2] + falj*ntj[2];

          gb_force[0] *= 2.0*epsilon_sigmasq_r2inv;
          gb_force[1] *= 2.0*epsilon_sigmasq_r2inv;
          gb_force[2] *= 2.0*epsilon_sigmasq_r2inv;

          eng = epsilon_sigmasq_r2inv*uij;

          energy[ET_BEND] += eng;
          if (eflag) one_eng += eng;

          // Adding all forces
          fforce[0] += fpair*r12[0] + gb_force[0];
          fforce[1] += fpair*r12[1] + gb_force[1];
          fforce[2] += fpair*r12[2] + gb_force[2];

          // Calculate torque on particle i

          // Torque from derivative over (ru, normi)
          tor_thi = 2.0*epsilon_sigmasq_r2inv*thetai_gamma_r;
          MathExtra::cross3(ru12, normi, torque);
          ttor[0] += tor_thi*torque[0];
          ttor[1] += tor_thi*torque[1];
          ttor[2] += tor_thi*torque[2];

          // Torque from derivative over (ru, nti)
          tor_ali = -2.0*epsilon_sigmasq_r2inv*fali*r;
          MathExtra::cross3(ru12, nti, torque);
          ttor[0] += tor_ali*torque[0];
          ttor[1] += tor_ali*torque[1];
          ttor[2] += tor_ali*torque[2];

          // Calculate torque on particle j

          // Torque from derivative over (ru, normj)
          tor_thj = 2.0*epsilon_sigmasq_r2inv*thetaj_gamma_r;
          MathExtra::cross3(ru12, normj, torque);
          rtor[0] += tor_thj*torque[0];
          rtor[1] += tor_thj*torque[1];
          rtor[2] += tor_thj*torque[2];

          // Torque from derivative over (ru, ntj)
          tor_alj = -2.0*epsilon_sigmasq_r2inv*falj*r;
          MathExtra::cross3(ru12, ntj, torque);
          rtor[0] += tor_alj*torque[0];
          rtor[1] += tor_alj*torque[1];
          rtor[2] += tor_alj*torque[2];

          // SPAM variable dynamics
          // Calculating -dphi/dt for phi_b only
          if (spam_flag && b_spam_force) {

            dpb = 0.0;
            // Contribution from the derivative of prefactor
            if (b_K_force) dpb += K_force_factor*bend_epsilon[itype][jtype]*bend_k0[itype][jtype]*lipid_factor*r2inv*uij;
            // Contribution from the derivative of gamma
            if (b_fij_force) dpb -= fij_force_factor*epsilon_sigmasq_r2inv*bend_gamma_epsilon[itype][jtype]*(thetai_gamma_r - thetaj_gamma_r)*alpha_sq*r;
            
            if (b_spam_force_i) idphi[1] += dpb;
            if (b_spam_force_j) jdphi[1] += dpb;
          }

          // Calculate -dphi/dt for phi_m only
          // Contribution from the derivative of prefactor
          if (spam_flag && m_spam_force) {
            dpm = -0.25*(lipid_lambda[itype][jtype] - 1.0)*bend_epsilon[itype][jtype]*(1 - bend_k0[itype][jtype]*K_factor)*r2inv*uij;

            if (iphi_m<=1.0 && iphi_m>=-1.0) idphi[0] += dpb;
            if (jphi_m<=1.0 && jphi_m>=-1.0) jdphi[0] += dpb;
          }
        }


        // Oligomerization potential
        // Torques and energy are calculated
        if (spam_flag && olig_pot_flag[itype][jtype]) {
          if (iphi_b<-0.8 && jphi_b<-0.8) {
            epsilon = olig_epsilon[itype][jtype];
            if (aolig[itype][jtype]==2) epsilon *= MathExtra::dot3(nti, ntj);
//            dot = MathExtra::dot3(nti, ntj);
//            epsilon *= pow(dot, aolig[itype][jtype]-1.0)
            MathExtra::cross3(nti, ntj, torque);
            torque[0] *= epsilon;
            torque[1] *= epsilon;
            torque[2] *= epsilon;

            eng = -epsilon*MathExtra::dot3(nti, ntj)/((double)aolig[itype][jtype]);

            energy[ET_OLIG] += eng;
            if (eflag) one_eng += eng;

            ttor[0] += torque[0];
            ttor[1] += torque[1];
            ttor[2] += torque[2];

            rtor[0] -= torque[0];
            rtor[1] -= torque[1];
            rtor[2] -= torque[2];
          }
        }

        // Two-body term of Intrinsic Curvature Coupling potential
        // Force, torque, dphi_b and energy are calculated
        // The one-body term is treated in the single loop bellow
        if (spam_flag && ic_pot_flag[itype] && ic_pot_flag[jtype]) {

          epsilon = 0.0;
          if (n_c[i]>0.0) epsilon += ic_lambda_k[itype]*MIN(MAX(iphi_b,-1.0),1.0)/n_c[i];
          if (n_c[j]>0.0) epsilon += ic_lambda_k[jtype]*MIN(MAX(jphi_b,-1.0),1.0)/n_c[j];

//          epsilon = ic_lambda_k[itype]*MIN(MAX(iphi_b,-1.0),1.0) + ic_lambda_k[jtype]*MIN(MAX(jphi_b,-1.0),1.0);

          rinv = 1.0/r;
          r2inv = 1.0/rsq;

          ru12[0] = r12[0]*rinv;
          ru12[1] = r12[1]*rinv;
          ru12[2] = r12[2]*rinv;

          thetai = MathExtra::dot3(normi,ru12);
          thetaj = MathExtra::dot3(normj,ru12);

          alphai = MathExtra::dot3(nti,ru12);
          alphaj = MathExtra::dot3(ntj,ru12);

          alpha_sq = alphai*alphai + alphaj*alphaj;

          // gamma_epsilon = (1/8)*(gamma_epsilon_i*fi + gamma_epsilon_j*fj)
          if (Gpd_flag) {
            fi = 0.5*(MAX(MIN(-iphi_b,1.0),-1.0) + 1.0);
            fj = 0.5*(MAX(MIN(-jphi_b,1.0),-1.0) + 1.0);
          } else {
            fi = MAX(MIN(-iphi_b,1.0),0.0);
            fj = MAX(MIN(-jphi_b,1.0),0.0);
          }
          gamma_epsilon = 0.125*(ic_gamma_epsilon[itype]*fi + ic_gamma_epsilon[jtype]*fj);

          gamma = gamma_epsilon*alpha_sq;
          gamma_r = gamma*r;

          thetai_gamma_r = thetai - gamma_r;
          thetaj_gamma_r = thetaj + gamma_r;

          uij = thetai_gamma_r*thetai_gamma_r + thetaj_gamma_r*thetaj_gamma_r;

          // Spherical force
          fpair = -(uij + thetai_gamma_r*thetai + thetaj_gamma_r*thetaj)*r2inv;
          fpair -= (thetai_gamma_r - thetaj_gamma_r)*(gamma - 2.0*gamma_epsilon*alpha_sq)*rinv;
          fpair *= 2.0*epsilon*r2inv;

          fthi = thetai_gamma_r*rinv;
          fthj = thetaj_gamma_r*rinv;

          fali = -2.0*gamma_epsilon*(thetai_gamma_r - thetaj_gamma_r)*alphai;
          falj = -2.0*gamma_epsilon*(thetai_gamma_r - thetaj_gamma_r)*alphaj;

          // Aspherical forces
          gb_force[0] = fthi*normi[0] + fthj*normj[0] + fali*nti[0] + falj*ntj[0];
          gb_force[1] = fthi*normi[1] + fthj*normj[1] + fali*nti[1] + falj*ntj[1];
          gb_force[2] = fthi*normi[2] + fthj*normj[2] + fali*nti[2] + falj*ntj[2];

          gb_force[0] *= 2.0*epsilon*r2inv;
          gb_force[1] *= 2.0*epsilon*r2inv;
          gb_force[2] *= 2.0*epsilon*r2inv;

          eng = -epsilon*r2inv*uij;

          energy[ET_IC] += eng;
          if (eflag) one_eng += eng;

          // Adding all forces
          fforce[0] += fpair*r12[0] + gb_force[0];
          fforce[1] += fpair*r12[1] + gb_force[1];
          fforce[2] += fpair*r12[2] + gb_force[2];

          // Calculate torque on particle i

          // Torque from derivative over (ru, normi)
          tor_thi = -2.0*epsilon*r2inv*thetai_gamma_r;
          MathExtra::cross3(ru12, normi, torque);
          ttor[0] += tor_thi*torque[0];
          ttor[1] += tor_thi*torque[1];
          ttor[2] += tor_thi*torque[2];

          // Torque from derivative over (ru, nti)
          tor_ali = -2.0*epsilon*r2inv*fali*r;
          MathExtra::cross3(ru12, nti, torque);
          ttor[0] += tor_ali*torque[0];
          ttor[1] += tor_ali*torque[1];
          ttor[2] += tor_ali*torque[2];

          // Calculate torque on particle j

          // Torque from derivative over (ru, normj)
          tor_thj = -2.0*epsilon*r2inv*thetaj_gamma_r;
          MathExtra::cross3(ru12, normj, torque);
          rtor[0] += tor_thj*torque[0];
          rtor[1] += tor_thj*torque[1];
          rtor[2] += tor_thj*torque[2];

          // Torque from derivative over (ru, ntj)
          tor_alj = -2.0*epsilon*r2inv*falj*r;
          MathExtra::cross3(ru12, ntj, torque);
          rtor[0] += tor_alj*torque[0];
          rtor[1] += tor_alj*torque[1];
          rtor[2] += tor_alj*torque[2];

          // SPAM variable dynamics
          // Calculating -dphi/dt for phi_b only
          if (b_spam_force2) {

            // There are two contributions:
            // from the derivative of lambda*phi in front,
            // and from the derivative of gamma
            dpb = 0.25*epsilon*r2inv*(thetai_gamma_r - thetaj_gamma_r)*alpha_sq*r;
            if (iphi_b<=1.0 && iphi_b>=-1.0 && n_c[i]>0.0) {

              idphi[1] += ic_lambda_k[itype]*uij*r2inv/n_c[i];
              if (Gpd_flag) idphi[1] += 0.5*ic_gamma_epsilon[itype]*dpb;
              else if (iphi_b<=0.0) idphi[1] += ic_gamma_epsilon[itype]*dpb;
            }
            if (newton_pair || j < nlocal) {
              if (jphi_b<=1.0 && jphi_b>=-1.0 && n_c[j]>0.0) {
                jdphi[1] += ic_lambda_k[jtype]*uij*r2inv/n_c[j];
                if (Gpd_flag) jdphi[1] += 0.5*ic_gamma_epsilon[jtype]*dpb;
                else if (jphi_b<=0.0) jdphi[1] += ic_gamma_epsilon[jtype]*dpb;
              }
            }
          }
        }

        // Membrain Composition potential
        // Contribution from gradient of phi_m term
        // Only dphi_m is caclulated
        // Energy and one-body well potential are calculated bellow
        // Unlike the older EM2 code, rho_m (rather then rho_b) is used for calculations here
        if (spam_flag && mem_comp_pot_flag[itype] && mem_comp_pot_flag[jtype] && rsq<comp_rcutsq) {
          sigma = r*comp_rcut_inv;
          rho_ave_inv = 1.0/(0.5*(rho_m[i]+rho_m[j]));
          dw = -12.0*(1.0 - sigma)*(1.0 - sigma)*comp_rcut_inv*comp_rcut_inv;
          dpm = 2.0*rho_ave_inv*dw;
          // dot(xi_m_i+xi_m_j, r12)
          dpm *= (xi_m[i][0]+xi_m[j][0])*r12[0] + (xi_m[i][1]+xi_m[j][1])*r12[1] + (xi_m[i][2]+xi_m[j][2])*r12[2];

          idphi[0] += mem_comp_xi_epsilon[itype]*mass[jtype]*dpm;
          if (newton_pair || j < nlocal) {
            jdphi[0] -= mem_comp_xi_epsilon[jtype]*mass[itype]*dpm;
          }
        }

        // Protein Composition potential
        // Contribution from gradient of phi_b term
        // Only dphi_b is caclulated
        // Energy and one-body well potential are calculated bellow
        if (spam_flag && prot_comp_pot_flag[itype] && prot_comp_pot_flag[jtype] && rsq<comp_rcutsq) {
          sigma = r*comp_rcut_inv;
          rho_ave_inv = 1.0/(0.5*(rho_b[i]+rho_b[j]));
          dw = -12.0*(1.0 - sigma)*(1.0 - sigma)*comp_rcut_inv*comp_rcut_inv;
          dpb = 2.0*rho_ave_inv*dw;
          // dot(xi_b_i+xi_b_j, r12)
          dpb *= (xi_b[i][0]+xi_b[j][0])*r12[0] + (xi_b[i][1]+xi_b[j][1])*r12[1] + (xi_b[i][2]+xi_b[j][2])*r12[2];

          idphi[1] += prot_comp_xi_epsilon[itype]*mass[jtype]*dpb;
          if (newton_pair || j < nlocal) {
            jdphi[1] -= prot_comp_xi_epsilon[jtype]*mass[itype]*dpb;
          }
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
        tor[i][1] += ttor[1];
        tor[i][2] += ttor[2];
        dphi[i][0] += spam_gamma*idphi[0];
        dphi[i][1] += spam_gamma*idphi[1];

        if (newton_pair || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
          f[j][1] -= fforce[1];
          f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
          tor[j][1] += rtor[1];
          tor[j][2] += rtor[2];
          dphi[j][0] += spam_gamma*jdphi[0];
          dphi[j][1] += spam_gamma*jdphi[1];
        }

        if (eflag) evdwl = factor_lj*one_eng;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fforce[0],fforce[1],fforce[2],
                                 r12[0],r12[1],r12[2]);
      }
    }
  }

  sum_phi_m = sum_phi_b = 0.0;
  n_phi_m = n_phi_b = 0.0;
  pcov_sum_phi_b = 0.0;
  pcov_n_phi_b = 0.0;
  pcov_n_sol_phi_b = 0.0;

  // One-body terms
  // idphi here is multipied by spam_gamma in the end
  for (i=0; i < nlocal; i++) {
    itype = type[i];

    one_eng = 0.0;
    idphi[0] = idphi[1] = 0.0;

    // SPAM variables, phi array
    // 0 - membrain, 1 - protein
    if (spam_flag) {
      iphi_m = phi[i][0];
      iphi_b = phi[i][1];
    }
    
    // One-body term of Intrinsic Curvature Coupling potential
    // Only dphi_b and energy are calculated
    // The pair term is treated in the double loop above
    if (spam_flag && ic_pot_flag[itype]) {
      eng = ic_lambda_m[itype]*MIN(MAX(iphi_b,-1.0),1.0);
      energy[ET_IC] += eng;
      if (eflag) one_eng += eng;

      if (iphi_b>=-1.0 && iphi_b<=1.0) {
        dpb = -ic_lambda_m[itype];
        idphi[1] += dpb;
      }
    }

    // Composition Coupling potential
    // Only dphi_m, dphi_b, and energy are callculated
    if (spam_flag && cc_pot_flag[itype]) {
      zeta_m = 0.5*(iphi_m - 1.0) - cc_zeta0[itype];
      zeta_b = 0.5*(iphi_b - 1.0) - cc_zeta0[itype];
      dpm = 0.5*cc_epsilon[itype]*zeta_b;
      dpb = 0.5*cc_epsilon[itype]*zeta_m;

      eng = -cc_epsilon[itype]*zeta_m*zeta_b;

      energy[ET_CC] += eng;
      if (eflag) one_eng += eng;

      idphi[0] += dpm;
      idphi[1] += dpb;
    }

    // Membrain composition potential
    // phi_m^10/10 or general polynomial potential
    // dphi_m from one-body well potential and total energy are calculated
    // dphi_m from gradient is calculated in the double loop above
    if (spam_flag && mem_comp_pot_flag[itype]) {
      // phi_m gradient term energy
      eng = mem_comp_xi_epsilon[itype]*MathExtra::dot3(xi_m[i],xi_m[i]);

      if (mem_comp_npoly[itype]==0) { 
        // Standard phi_m^10/10 single-well potential

        dpm = mem_comp_epsilon[itype]*pow(iphi_m,9.0);
        eng += mem_comp_epsilon[itype]*0.1*pow(iphi_m,10.0);
      } else {
        // General polynomial potential

        dpm = 0.0;
        for (j=0;j<mem_comp_npoly[itype];j++) {
          dpm += mem_comp_poly_exp[itype][j]*mem_comp_poly_coeff[itype][j]*pow(iphi_m, mem_comp_poly_exp[itype][j]-1); 
          eng += mem_comp_poly_coeff[itype][j]*pow(iphi_m, mem_comp_poly_exp[itype][j]);
        }
      }

      energy[ET_MEM_COMP] += eng;
      if (eflag) one_eng += eng;

      idphi[0] -= dpm;
    }

    // Protein composition potential
    // phi_b^6 + phi_b^2
    // dphi_b from one-body well potential and total energy are calculated
    // dphi_b from gradient is calculated in the double loop above
    if (spam_flag && prot_comp_pot_flag[itype]) {
      phi2 = iphi_b*iphi_b;
      phi4 = phi2*phi2;
      dpb = prot_comp_epsilon[itype]*iphi_b*(6.0*phi4 + 2.0);

      // phi_b gradient term energy
      eng = prot_comp_xi_epsilon[itype]*MathExtra::dot3(xi_b[i],xi_b[i]);
      // Well potential
      eng += prot_comp_epsilon[itype]*phi2*(phi4 + 1.0);

      energy[ET_PROT_COMP] += eng;
      if (eflag) one_eng += eng;

      idphi[1] -= dpb;
    }

    // Sum over phi_m in order to propagate
    // composition-stat for membrain composition
    if (spam_flag && mem_stat_flag[itype]) {
      sum_phi_m += iphi_m;
      n_phi_m += 1.0;
    }

    // Sum over phi_b in order to propagate
    // composition-stat for protein composition
    if (spam_flag && prot_stat_flag[itype]) {
      sum_phi_b += iphi_b;
      n_phi_b += 1.0;
    }

    // Sum over phi_b on the membrane
    // To be used for Protein Coverage potential
    if (spam_flag && pcov_pot_flag[itype]) {
      if (pcov_pot_flag[itype]==1) {
        pcov_sum_phi_b += MIN(MAX(iphi_b,-1.0),1.0);
        pcov_n_phi_b += 1.0;
      } else {
        pcov_n_sol_phi_b += 1.0;
      }
    }

    dphi[i][0] += spam_gamma*idphi[0];
    dphi[i][1] += spam_gamma*idphi[1];

    if (eflag) evdwl = factor_lj*one_eng;

    if (eflag_global) eng_vdwl += evdwl;
    if (eflag_atom) eatom[i] += evdwl;
  }

  MPI_Allreduce(&sum_phi_m,&sum_phi_m_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&sum_phi_b,&sum_phi_b_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&n_phi_m,&n_phi_m_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&n_phi_b,&n_phi_b_all,1,MPI_DOUBLE,MPI_SUM,world);

  MPI_Allreduce(&pcov_sum_phi_b,&pcov_sum_phi_b_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&pcov_n_phi_b,&pcov_n_phi_b_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&pcov_n_sol_phi_b,&pcov_n_sol_phi_b_all,1,MPI_DOUBLE,MPI_SUM,world);

  // Protein Coverage potential
  // Energy calculation
  if (spam_flag && pcov_n_phi_b_all>0.0) {
//    if (comm->me==0) printf("step: %d pcov_sum_phi_b_all %f pcov_n_phi_b_all %f\n", update->ntimestep, pcov_sum_phi_b_all, pcov_n_phi_b_all);

    pcov_sum_phi_b_all *= 0.5/pcov_n_phi_b_all;

    if (comm->me==0) {
      eng = pcov_epsilon*pow(pcov_sum_phi_b_all - 0.5 + pcov_eta0, 2.0);
      energy[ET_PCOV] += eng;
      
//      printf("step: %d pcov_sum_phi_b_all %f pcov_n_phi_b_all %f eng %f\n", update->ntimestep, pcov_sum_phi_b_all, pcov_n_phi_b_all, eng);
    
      if (eflag) evdwl = factor_lj*eng;

      if (eflag_global) eng_vdwl += evdwl;
      if (eflag_atom) eatom[i] += evdwl;
    }
  }

  // Propogate mem_stat and prot_stat
  // stat += epsilon*sum_phi/n_phi

  if (n_phi_m_all>0.0) mem_stat += update->dt*mem_stat_epsilon*sum_phi_m_all/n_phi_m_all;
  if (n_phi_b_all>0.0) prot_stat += update->dt*prot_stat_epsilon*sum_phi_b_all/n_phi_b_all;

  // Second one-body loop
  // idphi here is NOT multipied by spam_gamma
  for (i=0; i < nlocal; i++) {
    itype = type[i];

    idphi[0] = idphi[1] = 0.0;

    if (spam_flag) {
      iphi_m = phi[i][0];
      iphi_b = phi[i][1];

      // Subtract stat from dphi
      if (mem_stat_flag[itype]) idphi[0] -= mem_stat;
      if (prot_stat_flag[itype]) idphi[1] -= prot_stat;

      // Calculating Protein Coverage potential
      // Only dphi_b will be calculated
      if (pcov_pot_flag[itype] && pcov_n_phi_b_all>0.0) {
        if (iphi_b>-1.0 && iphi_b<1.0) {
          idphi[1] -= pcov_epsilon*(pcov_sum_phi_b_all - 0.5 + pcov_eta0)/pcov_n_phi_b_all;
//          if (atom->tag[i]==1) printf("step: %d dphi: %f\n",update->ntimestep,pcov_epsilon*(pcov_sum_phi_b_all - 0.5 + pcov_eta0)/pcov_n_phi_b_all);
          //idphi[1] -= spam_gamma*pcov_epsilon*(pcov_sum_phi_b_all - 0.5 + pcov_eta0)/pcov_n_phi_b_all;
        }
      }

      // Velocity times composition gradient
      // Possibly need to move to integration
      if (flow_term_flag) {
        if (mem_comp_pot_flag[itype]) idphi[0] += MathExtra::dot3(v[i],xi_m[i]);
        if (prot_comp_pot_flag[itype]) idphi[1] += MathExtra::dot3(v[i],xi_b[i]);
      }
    }

    dphi[i][0] += idphi[0];
    dphi[i][1] += idphi[1];
  }

  // Sum energies across terms and across processors

  for (int i=1;i<nEnergyTerms;++i) energy[ET_TOTAL] += energy[i];
  MPI_Allreduce(energy,energy_all,nEnergyTerms,MPI_DOUBLE,MPI_SUM,world);


  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEM2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(lj216_pot_flag,n+1,n+1,"pair:lj216_pot_flag");
  memory->create(lj216_epsilon,n+1,n+1,"pair:lj216_epsilon");
  memory->create(lj216_sigma,n+1,n+1,"pair:lj216_sigma");
  memory->create(lj216_k0,n+1,n+1,"pair:lj216_k0");
  memory->create(lj612_1,n+1,n+1,"pair:lj612_1");
  memory->create(lj612_2,n+1,n+1,"pair:lj612_2");
  memory->create(lj612_3,n+1,n+1,"pair:lj612_3");
  memory->create(lj612_4,n+1,n+1,"pair:lj612_4");
  memory->create(lj612_pot_flag,n+1,n+1,"pair:lj612_pot_flag");
  memory->create(lj612_epsilon,n+1,n+1,"pair:lj612_epsilon");
  memory->create(lj612_sigma,n+1,n+1,"pair:lj612_sigma");
  memory->create(lj612_k0,n+1,n+1,"pair:lj612_k0");
  memory->create(lucy_pot_flag,n+1,n+1,"pair:lucy_pot_flag");
  memory->create(lucy_epsilon,n+1,n+1,"pair:lucy_epsilon");
  memory->create(lucy_sigma,n+1,n+1,"pair:lucy_sigma");
  memory->create(lucy_sigma_inv,n+1,n+1,"pair:lucy_sigma_inv");
  memory->create(lucy_sigmasq_inv,n+1,n+1,"pair:lucy_sigmasq_inv");
  memory->create(lucy_table_pot_flag,n+1,n+1,"pair:lucy_table_pot_flag");
  memory->create(lucy_table_epsilon,n+1,n+1,"pair:lucy_table_epsilon");
  memory->create(lucy_table_sigma,n+1,n+1,"pair:lucy_table_sigma");
  memory->create(ex12_pot_flag,n+1,n+1,"pair:ex12_pot_flag");
  memory->create(ex12_epsilon,n+1,n+1,"pair:ex12_epsilon");
  memory->create(ex12_sigma,n+1,n+1,"pair:ex12_sigma");
  memory->create(ex12_epsig,n+1,n+1,"pair:ex12_epsig");
  memory->create(gauss_pot_flag,n+1,n+1,"pair:gauss_pot_flag");
  memory->create(gauss_epsilon,n+1,n+1,"pair:gauss_epsilon");
  memory->create(gauss_sigma,n+1,n+1,"pair:gauss_sigma");
  memory->create(gauss_sigmasq_inv,n+1,n+1,"pair:gauss_sigmasq_inv");
  memory->create(gauss_r0,n+1,n+1,"pair:gauss_r0");
  memory->create(bend_pot_flag,n+1,n+1,"pair:bend_pot_flag");
  memory->create(bend_epsilon,n+1,n+1,"pair:bend_epsilon");
  memory->create(bend_k0,n+1,n+1,"pair:bend_k0");
  memory->create(bend_gamma_epsilon,n+1,n+1,"pair:bend_gamma_epsilon");
  memory->create(lipid_factor_flag,n+1,n+1,"pair:lipid_factor_flag");
  memory->create(lipid_lambda,n+1,n+1,"pair:lipid_lambda");
  memory->create(olig_pot_flag,n+1,n+1,"pair:olig_pot_flag");
  memory->create(olig_epsilon,n+1,n+1,"pair:olig_epsilon");
  memory->create(aolig,n+1,n+1,"pair:aolig");
  memory->create(ic_pot_flag,n+1,"pair:ic_pot_flag");
  memory->create(ic_lambda_m,n+1,"pair:lambda_m");
  memory->create(ic_lambda_k,n+1,"pair:lambda_k");
  memory->create(ic_gamma_epsilon,n+1,"pair:ic_gamma_epsilon");
  memory->create(cc_pot_flag,n+1,"pair:cc_pot_flag");
  memory->create(cc_epsilon,n+1,"pair:cc_epsilon");
  memory->create(cc_zeta0,n+1,"pair:cc_zeta0");
  memory->create(mem_comp_pot_flag,n+1,"pair:mem_comp_pot_flag");
  memory->create(mem_comp_epsilon,n+1,"pair:mem_comp_epsilon");
  memory->create(mem_comp_xi_epsilon,n+1,"pair:mem_comp_xi_epsilon");
  memory->create(mem_comp_npoly,n+1,"pair:mem_comp_npoly");
//  memory->create(mem_comp_poly_exp,n+1,1,"pair:mem_comp_poly_exp");
//  memory->create(mem_comp_poly_coeff,n+1,1,"pair:mem_comp_poly_coeff");
  memory->create(prot_comp_pot_flag,n+1,"pair:prot_comp_pot_flag");
  memory->create(prot_comp_epsilon,n+1,"pair:prot_comp_epsilon");
  memory->create(prot_comp_xi_epsilon,n+1,"pair:prot_comp_xi_epsilon");
  memory->create(mem_stat_flag,n+1,"pair:mem_stat_flag");
  memory->create(prot_stat_flag,n+1,"pair:prot_stat_flag");
  memory->create(pcov_pot_flag,n+1,"pair:pcov_pot_flag");

  for (int i = 1; i <= n; i++) {
    ic_pot_flag[i] = 0;
    cc_pot_flag[i] = 0;
    mem_comp_pot_flag[i] = 0;
    mem_comp_npoly[i] = 0;
    prot_comp_pot_flag[i] = 0;
    mem_stat_flag[i] = 0;
    prot_stat_flag[i] = 0;
    pcov_pot_flag[i] = 0;
    for (int j = 1; j <= n; j++) {
      lj216_pot_flag[i][j] = 0;
      lj612_pot_flag[i][j] = 0;
      lucy_pot_flag[i][j] = 0;
      lucy_table_pot_flag[i][j] = 0;
      ex12_pot_flag[i][j] = 0;
      gauss_pot_flag[i][j] = 0;
      bend_pot_flag[i][j] = 0;
      lipid_factor_flag[i][j] = 0;
      olig_pot_flag[i][j] = 0;
    }
  }

  mem_comp_poly_exp = new double*[n+1];
  mem_comp_poly_coeff = new double*[n+1];
  for (int i = 1; i <= n; i++) {
    mem_comp_poly_exp[i] = NULL;
    mem_comp_poly_coeff[i] = NULL;
  }

  lucy_etable1 = (double ***)malloc((n+1)*sizeof(double**));
  lucy_etable2 = (double ***)malloc((n+1)*sizeof(double**));
  lucy_ftable1 = (double ***)malloc((n+1)*sizeof(double**));
  lucy_ftable2 = (double ***)malloc((n+1)*sizeof(double**));
  for (int i = 0; i <= n; i++) {
    lucy_etable1[i] = (double **) malloc((n+1)*sizeof(double *));
    lucy_etable2[i] = (double **) malloc((n+1)*sizeof(double *));
    lucy_ftable1[i] = (double **) malloc((n+1)*sizeof(double *));
    lucy_ftable2[i] = (double **) malloc((n+1)*sizeof(double *));
    for (int j = 1; j <= n; j++) {
      lucy_etable1[i][j] = NULL;
      lucy_etable2[i][j] = NULL;
      lucy_ftable1[i][j] = NULL;
      lucy_ftable2[i][j] = NULL;
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEM2::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  spam_flag = 0; 
  flow_term_flag = 0;

  comp_flag = 0;
  mem_stat = 0.0;
  prot_stat = 0.0;

  // Flags for protein concentration dependence of Ks and gamma.
  Kpd_flag = 0;
  Gpd_flag = 0;

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  read_par_flag = true;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEM2::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int len = strlen(arg[2]) + 1;
  parfile = new char[len];
  strcpy(parfile, arg[2]);

  double cut_one = cut_global;
  if (narg >=4) cut_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  if (read_par_flag) read_parameters();
  compute_lucy_potential_table();
}

/* ----------------------------------------------------------------------*/

void PairEM2::read_parameters()
{
  read_par_flag = false;

  int i, j;
  int ilo,ihi,jlo,jhi;
  double epsilon_val, epsilon2_val, sigma_val, k0_val, geps_val;
  double lambda_m_val, lambda_k_val, zeta0_val, r0_val;
  int aolig_val;
  int flag1, flag2;
  int npol_val, npol;
  double coeff_one;
  double *pol_exp=NULL;
  double *pol_coeff=NULL;

  FILE *file;
  int file_state, narg=0, iarg_line=-1;
  char ln[1024], *line, *arg[16];
  enum File_States{FS_NONE=0, FS_LJ216, FS_LJ612, FS_LUCY, FS_LUCY_TB, FS_EX12, FS_GAUSS, FS_BEND, FS_OLIG, FS_COUP_IC, FS_COUP_CC, FS_MEM_COMP, FS_MEM_COMP_POLY, FS_BAR_COMP, FS_COMP_CUT, FS_COMP_STAT, FS_FLOW, FS_SPAM, FS_PCOV, FS_LIP_FACTOR, FS_K_FLAGS};

  int me = comm->me;

  // Open parameter file and read it line by line

  file = fopen(parfile,"r");
  if (!file) error->all(FLERR,"Pair_style EM2: Error opening coefficient file!");
  file_state = FS_NONE;
  while ( fgets ( ln, sizeof ln, file ) != NULL ) {
    line = trim(ln);

    if (line[0]=='#') continue;
    if (line[0]=='[') file_state = FS_NONE;
    if (isEmptyString(line)) { file_state = FS_NONE; continue; }

    // Split the parameter line into strings
    if (file_state!=FS_NONE) {
      narg = 0;
      arg[narg] = strtok(line, " \t\n");
      while ( arg[narg]!=NULL ) {
        narg++;
        if (narg>=16) error->all(FLERR,"Pair_style EM2: Parameter line is too long!");
        arg[narg] = strtok(NULL, " \t\n");
      }
    }

    iarg_line++;

    switch (file_state) {
    case FS_LJ216:
      if (narg!=5) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (LJ216)");
      if (me==0) print_log("LJ_2-16 potential flag on\n");
      epsilon_val = atof(arg[2]);
      k0_val = atof(arg[3]);
      sigma_val = atof(arg[4]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          lj216_pot_flag[i][j] = 1;
          lj216_epsilon[i][j] = epsilon_val;
          lj216_k0[i][j] = k0_val;
          lj216_sigma[i][j] = sigma_val;
        }
      }
      break;
    case FS_LJ612:
      if (narg!=5) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (LJ612)");
      if (me==0) print_log("LJ_2-16 potential flag on\n");
      epsilon_val = atof(arg[2]);
      k0_val = atof(arg[3]);
      sigma_val = atof(arg[4]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          lj612_pot_flag[i][j] = 1;
          lj612_epsilon[i][j] = epsilon_val;
          lj612_k0[i][j] = k0_val;
          lj612_sigma[i][j] = sigma_val;
        }
      }
      break;
    case FS_LUCY:
      if (narg!=4) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Lucy)");
      if (me==0) print_log("Lucy excluded volume potential flag on\n");
      epsilon_val = atof(arg[2]);
      sigma_val = atof(arg[3]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          lucy_pot_flag[i][j] = 1;
          lucy_epsilon[i][j] = epsilon_val;
          lucy_sigma[i][j] = sigma_val;
        }
      }
      break;
    case FS_LUCY_TB:
      if (iarg_line==0) {
        if (narg!=1) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Lucy Table Header)");
        if (me==0) print_log("Reading Lucy Table Header\n");
        lucy_tb_dr = atof(arg[0]);
        lucy_tb_dr_inv = 1.0/lucy_tb_dr;
      } else {
        if (narg!=4) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Lucy Table)");
        if (me==0) print_log("Lucy Table excluded volume potential flag on\n");
        epsilon_val = atof(arg[2]);
        sigma_val = atof(arg[3]);

        force->bounds(arg[0],atom->ntypes,ilo,ihi);
        force->bounds(arg[1],atom->ntypes,jlo,jhi);

        for (i = ilo; i <= ihi; i++) {
          for (j = MAX(jlo,i); j <= jhi; j++) {
            lucy_table_pot_flag[i][j] = 1;
            lucy_table_epsilon[i][j] = epsilon_val;
            lucy_table_sigma[i][j] = sigma_val;
          }
        }
      }
      break;
    case FS_EX12:
      if (narg!=4) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Excluded_12)");
      if (me==0) print_log("Excluded_12 potential flag on\n");
      epsilon_val = atof(arg[2]);
      sigma_val = atof(arg[3]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          ex12_pot_flag[i][j] = 1;
          ex12_epsilon[i][j] = epsilon_val;
          ex12_sigma[i][j] = sigma_val;
        }
      }
      break;
    case FS_GAUSS:
      if (narg!=5) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Gauss)");
      if (me==0) print_log("Gaussian potential flag on\n");
      epsilon_val = atof(arg[2]);
      sigma_val = atof(arg[3]);
      r0_val = atof(arg[4]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          gauss_pot_flag[i][j] = 1;
          gauss_epsilon[i][j] = epsilon_val;
          gauss_sigma[i][j] = sigma_val;
          gauss_r0[i][j] = r0_val;
        }
      }
      break;
    case FS_BEND:
      if (narg!=6) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Bending Pot)");
      if (me==0) print_log("Bending potential flag on\n");
      epsilon_val = atof(arg[2]);
      k0_val = atof(arg[3]);
      sigma_val = atof(arg[4]);
      geps_val = atof(arg[5]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          bend_pot_flag[i][j] = 1;
          bend_epsilon[i][j] = 4.0*epsilon_val*sigma_val*sigma_val; // 4*epsilon*sigma^2
          bend_k0[i][j] = k0_val;
          bend_gamma_epsilon[i][j] = geps_val;
        }
      }
      break;
    case FS_LIP_FACTOR:
      if (narg!=3) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Lipid Bending Rigidity Factor)");
      if (me==0) print_log("Lipid Bending Rigidity Factor flag on\n");
      epsilon_val = atof(arg[2]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          lipid_factor_flag[i][j] = 1;
          lipid_lambda[i][j] = epsilon_val;
        }
      }
      break;
    case FS_OLIG:
      if (narg!=4) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Oligomerization Pot)");
      if (me==0) print_log("Oligomerization potential flag on\n");
      epsilon_val = atof(arg[2]);
      aolig_val = atoi(arg[3]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);
      force->bounds(arg[1],atom->ntypes,jlo,jhi);

      for (i = ilo; i <= ihi; i++) {
        for (j = MAX(jlo,i); j <= jhi; j++) {
          olig_pot_flag[i][j] = 1;
          olig_epsilon[i][j] = epsilon_val;
          aolig[i][j] = aolig_val;
        }
      }
      break;
    case FS_COUP_IC:
      if (narg!=4) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Intrinsic Curvature Coupling Pot)");
      if (me==0) print_log("Intrinsic Curvature Coupling potential flag on\n");
      lambda_m_val = atof(arg[1]);
      lambda_k_val = atof(arg[2]);
      geps_val = atof(arg[3]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);

      for (i = ilo; i <= ihi; i++) {
        ic_pot_flag[i] = 1;
        ic_lambda_m[i] = lambda_m_val;
        ic_lambda_k[i] = lambda_k_val;
        ic_gamma_epsilon[i] = geps_val;
      }
      break;
    case FS_COUP_CC:
      if (narg!=3) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Composition Coupling Pot)");
      if (me==0) print_log("Composition Coupling potential flag on\n");
      epsilon_val = atof(arg[1]);
      zeta0_val = atof(arg[2]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);

      for (i = ilo; i <= ihi; i++) {
        cc_pot_flag[i] = 1;
        cc_epsilon[i] = epsilon_val;
        cc_zeta0[i] = zeta0_val;
      }
      break;
    case FS_MEM_COMP:
      if (narg!=3) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Membrane Composition Pot)");
      if (me==0) print_log("Membrane Composition potential flag on\n");
      epsilon_val = atof(arg[1]);
      epsilon2_val = atof(arg[2]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);

      for (i = ilo; i <= ihi; i++) {
        if (mem_comp_pot_flag[i]==1) error->all(FLERR,"Pair_style EM2: Repeated definition of parameters for a type (Membrane Composition Pot)");
        mem_comp_pot_flag[i] = 1;
        mem_comp_xi_epsilon[i] = epsilon_val;
        mem_comp_epsilon[i] = epsilon2_val;
      }
      break;
    case FS_MEM_COMP_POLY:
      if (narg<3) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Membrane Composition Poly Pot)");
      if (me==0) print_log("Membrane Composition Poly potential flag on\n");
      epsilon_val = atof(arg[1]);
      npol_val = atoi(arg[2]);
      if (npol_val<=0 || narg-3!=npol_val+1) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Membrane Composition Poly Pot)");

      memory->destroy(pol_exp);
      memory->destroy(pol_coeff);
      memory->create(pol_exp,npol_val+1,"pair:pol_exp");
      memory->create(pol_coeff,npol_val+1,"pair:pol_coeff");

      npol = 0;
      for (i=0;i<=npol_val;i++) {
        coeff_one = atof(arg[3+i]);
        if (coeff_one!=0.0) {
          pol_exp[npol] = (double)i;
          pol_coeff[npol] = coeff_one;
          npol++;
        }
      }
      if (npol<=0) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Membrane Composition Poly Pot)");

      force->bounds(arg[0],atom->ntypes,ilo,ihi);

      for (i = ilo; i <= ihi; i++) {
        if (mem_comp_pot_flag[i]==1) error->all(FLERR,"Pair_style EM2: Repeated definition of parameters for a type (Membrane Composition Poly Pot)");
        mem_comp_pot_flag[i] = 1;
        mem_comp_xi_epsilon[i] = epsilon_val;
        mem_comp_npoly[i] = npol;
        if (mem_comp_poly_exp[i]!=NULL) delete [] mem_comp_poly_exp[i];
        if (mem_comp_poly_coeff[i]!=NULL) delete [] mem_comp_poly_coeff[i];
        mem_comp_poly_exp[i] = new double[npol];
        mem_comp_poly_coeff[i] = new double[npol];
        for (j=0;j<npol;j++) {
          mem_comp_poly_exp[i][j] = pol_exp[j];
          mem_comp_poly_coeff[i][j] = pol_coeff[j];
        }
      }
      break;
    case FS_BAR_COMP:
      if (narg!=3) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Protein Composition Pot)");
      if (me==0) print_log("Protein Composition potential flag on\n");
      epsilon_val = atof(arg[1]);
      epsilon2_val = atof(arg[2]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);

      for (i = ilo; i <= ihi; i++) {
        prot_comp_pot_flag[i] = 1;
        prot_comp_xi_epsilon[i] = epsilon_val;
        prot_comp_epsilon[i] = epsilon2_val;
      }
      break;
    case FS_COMP_CUT:
      if (iarg_line!=0 || narg!=1) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Composition Cutoff)");
      if (me==0) print_log("Composition Cutoff flag on\n");
      comp_flag = 1;
      sigma_val = atof(arg[0]);
      comp_rcut = sigma_val;
      comp_rcutsq = sigma_val*sigma_val;
      comp_rcut_inv = 1.0/sigma_val;
      break;
    case FS_COMP_STAT:
      if (iarg_line==0) {
        if (narg!=2) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Composition Stat Header)");
        if (me==0) print_log("Reading Composition Stat flags\n");
        mem_stat_epsilon = atof(arg[0]);
        prot_stat_epsilon = atof(arg[1]);
      } else {
        if (narg!=3) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Composition Stat)");
        if (me==0) print_log("Reading Composition Stat flags\n");
        flag1 = atoi(arg[1]);
        flag2 = atoi(arg[2]);

        if ((flag1!=0 && flag1!=1) || (flag2!=0 && flag2!=1))
          error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Composition Stat)");

        force->bounds(arg[0],atom->ntypes,ilo,ihi);

        for (i = ilo; i <= ihi; i++) {
          mem_stat_flag[i] = flag1;
          prot_stat_flag[i] = flag2;
        }
      }
      break;
    case FS_PCOV:
      if (iarg_line==0) {
        if (narg!=2) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Protein Coverage Header)");
        if (me==0) print_log("Reading Protein Coverage Header\n");
        pcov_epsilon = atof(arg[0]);
        pcov_eta0 = atof(arg[1]);
      } else {
        if (narg!=2) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Protein Coverage)");
        if (me==0) print_log("Reading Protein Coverage parameters\n");
        flag1 = atoi(arg[1]);

        if (flag1!=0 && flag1!=1)
          error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Protein Coverage)");

        force->bounds(arg[0],atom->ntypes,ilo,ihi);

        for (i = ilo; i <= ihi; i++) {
          pcov_pot_flag[i] = flag1;
        }
      }
      break;
    case FS_FLOW:
      if (iarg_line!=0 || narg!=1) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Flow Term Flag)");
      if (me==0) print_log("Reading Flow Term flag\n");

      flag1 = atoi(arg[0]);
      if (flag1!=0 && flag1!=1) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Flow Term Flag)");

      flow_term_flag = flag1;
      break;
    case FS_SPAM:
      if (iarg_line!=0 || narg!=2) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (SPAM Flag)");
      if (me==0) print_log("Reading SPAM flag\n");

      flag1 = atoi(arg[0]);
      epsilon_val = atof(arg[1]);

      if (flag1!=0 && flag1!=1) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (SPAM Flag)");

      spam_flag = flag1;
      spam_gamma = epsilon_val;
      break;
    case FS_K_FLAGS:
      if (iarg_line!=0 || narg!=2) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (K-Flags)");
      if (me==0) print_log("Reading K-Flags\n");

      flag1 = atoi(arg[0]);
      flag2 = atoi(arg[1]);

      if ((flag1!=0 && flag1!=1) || (flag2!=0 && flag2!=1)) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (K-Flags)");

      Kpd_flag = flag1;
      Gpd_flag = flag2;
      break;
    case FS_NONE:
      iarg_line = -1;
      if (strcmp(line, "[Lennard-Jones_2-16]")==0)
        file_state = FS_LJ216;
      if (strcmp(line, "[Lennard-Jones_6-12]")==0)
        file_state = FS_LJ612;
      else if (strcmp(line, "[Lucy_Excluded_Volume]")==0)
        file_state = FS_LUCY;
      else if (strcmp(line, "[Lucy_Excluded_Volume_Table]")==0)
        file_state = FS_LUCY_TB;
      else if (strcmp(line, "[Excluded_12]")==0)
        file_state = FS_EX12;
      else if (strcmp(line, "[Gaussian]")==0)
        file_state = FS_GAUSS;
      else if (strcmp(line, "[EM2_Bending]")==0)
        file_state = FS_BEND;
      else if (strcmp(line, "[Lipid_Bending_Rigidity_Factor]")==0)
        file_state = FS_LIP_FACTOR;
      else if (strcmp(line, "[Oligomerization_Energy]")==0)
        file_state = FS_OLIG;
      else if (strcmp(line, "[Intrinsic_Curvature_Coupling]")==0)
        file_state = FS_COUP_IC;
      else if (strcmp(line, "[Composition_Coupling]")==0)
        file_state = FS_COUP_CC;
      else if (strcmp(line, "[Membrane_Composition]")==0)
        file_state = FS_MEM_COMP;
      else if (strcmp(line, "[Membrane_Composition_Poly]")==0)
        file_state = FS_MEM_COMP_POLY;
    else if (strcmp(line, "[Protein_Composition]")==0)
        file_state = FS_BAR_COMP;
      else if (strcmp(line, "[Composition_Cutoff]")==0)
        file_state = FS_COMP_CUT;
      else if (strcmp(line, "[Composition_Stat]")==0)
        file_state = FS_COMP_STAT;
      else if (strcmp(line, "[Protein_Coverage]")==0)
        file_state = FS_PCOV;
      else if (strcmp(line, "[Flow_Term_Flag]")==0)
        file_state = FS_FLOW;
      else if (strcmp(line, "[SPAM_Flag]")==0)
        file_state = FS_SPAM;
      else if (strcmp(line, "[K-Flags]")==0)
        file_state = FS_K_FLAGS;
      break;
    }
  }
  fclose(file);
  
  memory->destroy(pol_exp);
  memory->destroy(pol_coeff);
}

inline void PairEM2::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

char *PairEM2::ltrim(char *s)
{     
  while(isspace(*s)) s++;     
  return s; 
}  

char *PairEM2::rtrim(char *s)
{
  char* back;
  int len = strlen(s);

  if(len == 0)
    return(s); 

  back = s + len;     
  while(isspace(*--back));     
  *(back+1) = '\0';     
  return s; 
}  

char *PairEM2::trim(char *s)
{     
  return rtrim(ltrim(s));  
}

bool PairEM2::isEmptyString(char *str)
{
  int len = strlen(str);
  
  if (len==0) return true;
  
  for (int i=0;i<len;++i) {
    if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') return false;
  }

  return true;
}

/* ----------------------------------------------------------------------
 * Compute table for Lucy Potential
 * ---------------------------------------------------------------------- */
void PairEM2::compute_lucy_potential_table()
{
  int i, j, k, nbins;
  double r1, r2, v1, v2, f1, f2, sigma0, sigma, epsq2, epsilon;
  int n = atom->ntypes;

  for (i=1;i<=n;i++) {
    for (j=1;j<=n;j++) {
      if (lucy_etable1[i][j]) free(lucy_etable1[i][j]);
      if (lucy_etable2[i][j]) free(lucy_etable2[i][j]);
      if (lucy_ftable1[i][j]) free(lucy_ftable1[i][j]);
      if (lucy_ftable2[i][j]) free(lucy_ftable2[i][j]);

      lucy_etable1[i][j] = NULL;
      lucy_etable2[i][j] = NULL;
      lucy_ftable1[i][j] = NULL;
      lucy_ftable2[i][j] = NULL;
    }
  }

  for (i=1;i<=n;i++) {
    for (j=i;j<=n;j++) {
      if (lucy_table_pot_flag[i][j]) {
        epsilon = lucy_table_epsilon[i][j];
        sigma0 = lucy_table_sigma[i][j];

        nbins = (int)(sigma0*lucy_tb_dr_inv)+1;
        lucy_etable1[i][j] = (double *)malloc(nbins*sizeof(double));
        lucy_etable2[i][j] = (double *)malloc(nbins*sizeof(double));
        lucy_ftable1[i][j] = (double *)malloc(nbins*sizeof(double));
        lucy_ftable2[i][j] = (double *)malloc(nbins*sizeof(double));
        
        for (k=0;k<nbins;k++) {
          r1 = (double)k*lucy_tb_dr;
          r2 = (double)(k+1)*lucy_tb_dr;

          sigma = r1/sigma0;
          epsq2 = epsilon*(1-sigma)*(1-sigma);
          v1 = epsq2*(1-sigma)*(1+3*sigma);
          f1 = 12*epsq2/(sigma0*sigma0);

          sigma = r2/sigma0;
          epsq2 = epsilon*(1-sigma)*(1-sigma);
          v2 = epsq2*(1-sigma)*(1+3*sigma);
          f2 = 12*epsq2/(sigma0*sigma0);

          lucy_etable1[i][j][k] = (v2-v1)/(r2-r1);
          lucy_etable2[i][j][k] = (v1*r2 - v2*r1)/(r2-r1);
          lucy_ftable1[i][j][k] = (f2-f1)/(r2-r1);
          lucy_ftable2[i][j][k] = (f1*r2 - f2*r1)/(r2-r1);
        }

        if (i!=j) {
          lucy_etable1[j][i] = lucy_etable1[i][j];
          lucy_etable2[j][i] = lucy_etable2[i][j];
          lucy_ftable1[j][i] = lucy_ftable1[i][j];
          lucy_ftable2[j][i] = lucy_ftable2[i][j];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEM2::init_style()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) error->all(FLERR,"Pair em2 requires atom style em2");

  neighbor->request(this,instance_me);

  norm0[0] = 0.0;
  norm0[1] = 0.0;
  norm0[2] = 1.0;

  nt0[0] = 0.0;
  nt0[1] = 1.0;
  nt0[2] = 0.0;

  for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;

  // Check that composition cutoff was set
  // if any of the corresponding potentioal is on

  int i;
  int flag = 0;
  for (i=1; i <= atom->ntypes; i++) {
    if (mem_comp_pot_flag[i] || prot_comp_pot_flag[i]) {
      flag = 1;
      break;
    }
  }

  lucy_tb_dr_inv = 1.0/lucy_tb_dr;

  if (flag && !comp_flag) error->all(FLERR,"Pair em2 requires composition cutoff if any corresponding potentials are on");
  if (!flag) comp_flag = 0;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEM2::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  if (lj216_pot_flag[i][j]) {
    lj1[i][j] = 64.0 * lj216_epsilon[i][j] * pow(lj216_sigma[i][j],16.0);
    lj2[i][j] = 8.0 * lj216_epsilon[i][j] * pow(lj216_sigma[i][j],2.0);
    lj3[i][j] = 4.0 * lj216_epsilon[i][j] * pow(lj216_sigma[i][j],16.0);
    lj4[i][j] = 4.0 * lj216_epsilon[i][j] * pow(lj216_sigma[i][j],2.0);

    if (offset_flag) {
      double ratio = lj216_sigma[i][j] / cut[i][j];
      offset[i][j] = 4.0 * lj216_epsilon[i][j] * (pow(ratio,16.0) - pow(ratio,2.0));
    } else offset[i][j] = 0.0;

    lj216_pot_flag[j][i] = lj216_pot_flag[i][j];
    lj216_epsilon[j][i] = lj216_epsilon[i][j];
    lj216_sigma[j][i] = lj216_sigma[i][j];
    lj1[j][i] = lj1[i][j];
    lj2[j][i] = lj2[i][j];
    lj3[j][i] = lj3[i][j];
    lj4[j][i] = lj4[i][j];
    offset[j][i] = offset[i][j];
  }

  if (lj612_pot_flag[i][j]) {
    lj612_1[i][j] = 48.0 * lj612_epsilon[i][j] * pow(lj612_sigma[i][j],12.0);
    lj612_2[i][j] = 24.0 * lj612_epsilon[i][j] * pow(lj612_sigma[i][j],6.0);
    lj612_3[i][j] = 4.0 * lj612_epsilon[i][j] * pow(lj612_sigma[i][j],12.0);
    lj612_4[i][j] = 4.0 * lj612_epsilon[i][j] * pow(lj612_sigma[i][j],6.0);

    if (offset_flag) {
      double ratio = lj612_sigma[i][j] / cut[i][j];
      offset[i][j] = 4.0 * lj612_epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
    } else offset[i][j] = 0.0;

    lj612_pot_flag[j][i] = lj612_pot_flag[i][j];
    lj612_epsilon[j][i] = lj612_epsilon[i][j];
    lj612_sigma[j][i] = lj612_sigma[i][j];
    lj612_1[j][i] = lj612_1[i][j];
    lj612_2[j][i] = lj612_2[i][j];
    lj612_3[j][i] = lj612_3[i][j];
    lj612_4[j][i] = lj612_4[i][j];
    offset[j][i] = offset[i][j];
  }

  if (lucy_pot_flag[i][j]) {
    if (lucy_sigma[i][j]>cut[i][j]) error->all(FLERR,"Lucy Excluded Volume cutoff cannot be larger than overal cutoff for a pair type");

    lucy_sigma_inv[i][j] = 1.0/lucy_sigma[i][j];
    lucy_sigmasq_inv[i][j] = 1.0/pow(lucy_sigma[i][j], 2.0);

    lucy_pot_flag[j][i] = lucy_pot_flag[i][j];
    lucy_epsilon[j][i] = lucy_epsilon[i][j];
    lucy_sigma[j][i] = lucy_sigma[i][j];
    lucy_sigma_inv[j][i] = lucy_sigma_inv[i][j];
    lucy_sigmasq_inv[j][i] = lucy_sigmasq_inv[i][j];
  }

  if (lucy_table_pot_flag[i][j]) {
    if (lucy_table_sigma[i][j]>cut[i][j]) error->all(FLERR,"Lucy Table Excluded Volume cutoff cannot be larger than overal cutoff for a pair type");

    lucy_table_pot_flag[j][i] = lucy_table_pot_flag[i][j];
    lucy_table_epsilon[j][i] = lucy_table_epsilon[i][j];
    lucy_table_sigma[j][i] = lucy_table_sigma[i][j];
  }

  if (ex12_pot_flag[i][j]) {
    ex12_epsig[i][j] = ex12_epsilon[i][j]*pow(ex12_sigma[i][j],12);

    ex12_pot_flag[j][i] = ex12_pot_flag[i][j];
    ex12_epsilon[j][i] = ex12_epsilon[i][j];
    ex12_sigma[j][i] = ex12_sigma[i][j];
    ex12_epsig[j][i] = ex12_epsig[i][j];
  }

  if (gauss_pot_flag[i][j]) {
    gauss_sigmasq_inv[i][j] = 1.0/(gauss_sigma[i][j]*gauss_sigma[i][j]);

    gauss_pot_flag[j][i] = gauss_pot_flag[i][j];
    gauss_epsilon[j][i] = gauss_epsilon[i][j];
    gauss_sigma[j][i] = gauss_sigma[i][j];
    gauss_sigmasq_inv[j][i] = gauss_sigmasq_inv[i][j];
    gauss_r0[j][i] = gauss_r0[i][j];
  }

  if (bend_pot_flag[i][j]) {
    bend_pot_flag[j][i] = bend_pot_flag[i][j]; 
    bend_epsilon[j][i] = bend_epsilon[i][j];
    bend_k0[j][i] = bend_k0[i][j];
    bend_gamma_epsilon[j][i] = bend_gamma_epsilon[i][j];
  }

  if (lipid_factor_flag[i][j]) {
    lipid_factor_flag[j][i] = lipid_factor_flag[i][j];
    lipid_lambda[j][i] = lipid_lambda[i][j];
  }

  if (olig_pot_flag[i][j]) {
    if (aolig[i][j]!=1 && aolig[i][j]!=2) error->all(FLERR,"Oligomerization potential energy exponent should either be 1 or 2");

    olig_pot_flag[j][i] = olig_pot_flag[i][j];
    olig_epsilon[j][i] = olig_epsilon[i][j];
    aolig[j][i] = aolig[i][j];
  }

  if ( (mem_comp_pot_flag[i] && mem_comp_pot_flag[j]) || (prot_comp_pot_flag[i] && prot_comp_pot_flag[j]) )
    if (comp_rcut>cut[i][j]) error->all(FLERR,"Composition cutoff cannot be larger than overal cutoff for a pair");

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEM2::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;

  // Single loop
  for (i = 1; i <= atom->ntypes; i++) {
    fwrite(&ic_pot_flag[i],sizeof(int),1,fp);
    fwrite(&cc_pot_flag[i],sizeof(int),1,fp);
    fwrite(&mem_comp_pot_flag[i],sizeof(int),1,fp);
    fwrite(&mem_comp_npoly[i],sizeof(int),1,fp);
    fwrite(&prot_comp_pot_flag[i],sizeof(int),1,fp);
    fwrite(&mem_stat_flag[i],sizeof(int),1,fp);
    fwrite(&prot_stat_flag[i],sizeof(int),1,fp);
    fwrite(&pcov_pot_flag[i],sizeof(int),1,fp);
    if (ic_pot_flag[i]) {
      fwrite(&ic_lambda_m[i],sizeof(double),1,fp);
      fwrite(&ic_lambda_k[i],sizeof(double),1,fp);
      fwrite(&ic_gamma_epsilon[i],sizeof(double),1,fp);
    }
    if (cc_pot_flag[i]) {
      fwrite(&cc_epsilon[i],sizeof(double),1,fp);
      fwrite(&cc_zeta0[i],sizeof(double),1,fp);
    }
    if (mem_comp_pot_flag[i]) {
      fwrite(&mem_comp_epsilon[i],sizeof(double),1,fp);
      fwrite(&mem_comp_xi_epsilon[i],sizeof(double),1,fp);
      if (mem_comp_npoly[i]>0) {
        for (j=0;j<mem_comp_npoly[i];j++) {
          fwrite(&mem_comp_poly_exp[i][j],sizeof(double),1,fp);
          fwrite(&mem_comp_poly_coeff[i][j],sizeof(double),1,fp);
        }
      }
    }
    if (prot_comp_pot_flag[i]) {
      fwrite(&prot_comp_epsilon[i],sizeof(double),1,fp);
      fwrite(&prot_comp_xi_epsilon[i],sizeof(double),1,fp);
    }
  }

  // Double loop
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      fwrite(&lj216_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&lj612_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&lucy_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&lucy_table_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&ex12_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&gauss_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&bend_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&lipid_factor_flag[i][j],sizeof(int),1,fp);
      fwrite(&olig_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&cut[i][j],sizeof(double),1,fp);
      if (lj216_pot_flag[i][j]) {
        fwrite(&lj216_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&lj216_k0[i][j],sizeof(double),1,fp);
        fwrite(&lj216_sigma[i][j],sizeof(double),1,fp);
      }
      if (lj612_pot_flag[i][j]) {
        fwrite(&lj612_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&lj612_k0[i][j],sizeof(double),1,fp);
        fwrite(&lj612_sigma[i][j],sizeof(double),1,fp);
      }
      if (lucy_pot_flag[i][j]) {
        fwrite(&lucy_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&lucy_sigma[i][j],sizeof(double),1,fp);
      }
      if (lucy_table_pot_flag[i][j]) {
        fwrite(&lucy_table_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&lucy_table_sigma[i][j],sizeof(double),1,fp);
      }
      if (ex12_pot_flag[i][j]) {
        fwrite(&ex12_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&ex12_sigma[i][j],sizeof(double),1,fp);
      }
      if (gauss_pot_flag[i][j]) {
        fwrite(&gauss_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&gauss_sigma[i][j],sizeof(double),1,fp);
        fwrite(&gauss_r0[i][j],sizeof(double),1,fp);
      }
      if (bend_pot_flag[i][j]) {
        fwrite(&bend_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&bend_k0[i][j],sizeof(double),1,fp);
        fwrite(&bend_gamma_epsilon[i][j],sizeof(double),1,fp);
      }
      if (lipid_factor_flag[i][j]) { 
        fwrite(&lipid_lambda[i][j],sizeof(double),1,fp);
      }
      if (olig_pot_flag[i][j]) {
        fwrite(&olig_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&aolig[i][j],sizeof(int),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEM2::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  if (!allocated) allocate();

  int i,j;
  int me = comm->me;

  // Single loop
  for (i = 1; i <= atom->ntypes; i++) {
    if (me == 0) {
      fread(&ic_pot_flag[i],sizeof(int),1,fp);
      fread(&cc_pot_flag[i],sizeof(int),1,fp);
      fread(&mem_comp_pot_flag[i],sizeof(int),1,fp);
      fread(&mem_comp_npoly[i],sizeof(int),1,fp);
      fread(&prot_comp_pot_flag[i],sizeof(int),1,fp);
      fread(&mem_stat_flag[i],sizeof(int),1,fp);
      fread(&prot_stat_flag[i],sizeof(int),1,fp);
      fread(&pcov_pot_flag[i],sizeof(int),1,fp);
    }
    MPI_Bcast(&ic_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&cc_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&mem_comp_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&mem_comp_npoly[i],1,MPI_INT,0,world);
    MPI_Bcast(&prot_comp_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&mem_stat_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&prot_stat_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&pcov_pot_flag[i],1,MPI_INT,0,world);

    if (ic_pot_flag[i]) {
      if (me == 0) {
        fread(&ic_lambda_m[i],sizeof(double),1,fp);
        fread(&ic_lambda_k[i],sizeof(double),1,fp);
        fread(&ic_gamma_epsilon[i],sizeof(double),1,fp);
      }
      MPI_Bcast(&ic_lambda_m[i],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&ic_lambda_k[i],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&ic_gamma_epsilon[i],1,MPI_DOUBLE,0,world);
    }

    if (cc_pot_flag[i]) {
      if (me == 0) {
        fread(&cc_epsilon[i],sizeof(double),1,fp);
        fread(&cc_zeta0[i],sizeof(double),1,fp);
      }
      MPI_Bcast(&cc_epsilon[i],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&cc_zeta0[i],1,MPI_DOUBLE,0,world);
    }

    if (mem_comp_pot_flag[i]) {
      if (me == 0) {
        fread(&mem_comp_epsilon[i],sizeof(double),1,fp);
        fread(&mem_comp_xi_epsilon[i],sizeof(double),1,fp);
        if (mem_comp_npoly[i]>0) {
          for (j=0;j<mem_comp_npoly[i];j++) {
            fread(&mem_comp_poly_exp[i][j],sizeof(double),1,fp);
            fread(&mem_comp_poly_coeff[i][j],sizeof(double),1,fp);
          }
        }
      }
      MPI_Bcast(&mem_comp_epsilon[i],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&mem_comp_xi_epsilon[i],1,MPI_DOUBLE,0,world);
      if (mem_comp_npoly[i]>0) {
        for (j=0;j<mem_comp_npoly[i];j++) {
          MPI_Bcast(&mem_comp_poly_exp[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&mem_comp_poly_coeff[i][j],1,MPI_DOUBLE,0,world);
        }
      }
    }

    if (prot_comp_pot_flag[i]) {
      if (me == 0) {
        fread(&prot_comp_epsilon[i],sizeof(double),1,fp);
        fread(&prot_comp_xi_epsilon[i],sizeof(double),1,fp);
      }
      MPI_Bcast(&prot_comp_epsilon[i],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&prot_comp_xi_epsilon[i],1,MPI_DOUBLE,0,world);
    }
  }

  // Double loop
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);

      if (me == 0) {
        fread(&lj216_pot_flag[i][j],sizeof(int),1,fp);
        fread(&lj612_pot_flag[i][j],sizeof(int),1,fp);
        fread(&lucy_pot_flag[i][j],sizeof(int),1,fp);
        fread(&lucy_table_pot_flag[i][j],sizeof(int),1,fp);
        fread(&ex12_pot_flag[i][j],sizeof(int),1,fp);
        fread(&gauss_pot_flag[i][j],sizeof(int),1,fp);
        fread(&bend_pot_flag[i][j],sizeof(int),1,fp);
        fread(&lipid_factor_flag[i][j],sizeof(int),1,fp);
        fread(&olig_pot_flag[i][j],sizeof(int),1,fp);
      }
      MPI_Bcast(&lj216_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&lj612_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&lucy_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&lucy_table_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&ex12_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&gauss_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&bend_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&lipid_factor_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&olig_pot_flag[i][j],1,MPI_INT,0,world);

      if (me == 0) fread(&cut[i][j],sizeof(double),1,fp);
      MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);

      if (lj216_pot_flag[i][j]) {
        if (me == 0) {
          fread(&lj216_epsilon[i][j],sizeof(double),1,fp);
          fread(&lj216_k0[i][j],sizeof(double),1,fp);
          fread(&lj216_sigma[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&lj216_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lj216_k0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lj216_sigma[i][j],1,MPI_DOUBLE,0,world);
      }

      if (lj612_pot_flag[i][j]) {
        if (me == 0) {
          fread(&lj612_epsilon[i][j],sizeof(double),1,fp);
          fread(&lj612_k0[i][j],sizeof(double),1,fp);
          fread(&lj612_sigma[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&lj612_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lj612_k0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lj612_sigma[i][j],1,MPI_DOUBLE,0,world);
      }

      if (lucy_pot_flag[i][j]) {
        if (me == 0) {
          fread(&lucy_epsilon[i][j],sizeof(double),1,fp);
          fread(&lucy_sigma[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&lucy_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lucy_sigma[i][j],1,MPI_DOUBLE,0,world);
      }

      if (lucy_table_pot_flag[i][j]) {
        if (me == 0) {
          fread(&lucy_table_epsilon[i][j],sizeof(double),1,fp);
          fread(&lucy_table_sigma[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&lucy_table_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lucy_table_sigma[i][j],1,MPI_DOUBLE,0,world);
      }

      if (ex12_pot_flag[i][j]) {
        if (me == 0) {
          fread(&ex12_epsilon[i][j],sizeof(double),1,fp);
          fread(&ex12_sigma[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&ex12_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&ex12_sigma[i][j],1,MPI_DOUBLE,0,world);
      }

      if (gauss_pot_flag[i][j]) {
        if (me == 0) {
          fread(&gauss_epsilon[i][j],sizeof(double),1,fp);
          fread(&gauss_sigma[i][j],sizeof(double),1,fp);
          fread(&gauss_r0[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&gauss_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gauss_sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gauss_r0[i][j],1,MPI_DOUBLE,0,world);
      }

      if (bend_pot_flag[i][j]) {
        if (me == 0) {
          fread(&bend_epsilon[i][j],sizeof(double),1,fp);
          fread(&bend_k0[i][j],sizeof(double),1,fp);
          fread(&bend_gamma_epsilon[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&bend_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&bend_k0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&bend_gamma_epsilon[i][j],1,MPI_DOUBLE,0,world);
      }

      if (lipid_factor_flag[i][j]) {
        if (me == 0) {
          fread(&lipid_lambda[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&lipid_lambda[i][j],1,MPI_DOUBLE,0,world);
      }

      if (olig_pot_flag[i][j]) {
        if (me == 0) {
          fread(&olig_epsilon[i][j],sizeof(double),1,fp);
          fread(&aolig[i][j],sizeof(int),1,fp);
        }
        MPI_Bcast(&olig_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&aolig[i][j],1,MPI_INT,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEM2::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&spam_flag,sizeof(int),1,fp);
  fwrite(&flow_term_flag,sizeof(int),1,fp);
  fwrite(&Kpd_flag,sizeof(int),1,fp);
  fwrite(&Gpd_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&comp_flag,sizeof(int),1,fp);
  fwrite(&mem_stat_epsilon,sizeof(double),1,fp);
  fwrite(&prot_stat_epsilon,sizeof(double),1,fp);
  fwrite(&mem_stat,sizeof(double),1,fp);
  fwrite(&prot_stat,sizeof(double),1,fp);
  fwrite(&pcov_epsilon,sizeof(double),1,fp);
  fwrite(&pcov_eta0,sizeof(double),1,fp);
  fwrite(&spam_gamma,sizeof(double),1,fp);
  fwrite(&lucy_tb_dr,sizeof(double),1,fp);

  if (comp_flag) {
    fwrite(&comp_rcut,sizeof(double),1,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEM2::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&spam_flag,sizeof(int),1,fp);
    fread(&flow_term_flag,sizeof(int),1,fp);
    fread(&Kpd_flag,sizeof(int),1,fp);
    fread(&Gpd_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&comp_flag,sizeof(int),1,fp);
    fread(&mem_stat_epsilon,sizeof(double),1,fp);
    fread(&prot_stat_epsilon,sizeof(double),1,fp);
    fread(&mem_stat,sizeof(double),1,fp);
    fread(&prot_stat,sizeof(double),1,fp);
    fread(&pcov_epsilon,sizeof(double),1,fp);
    fread(&pcov_eta0,sizeof(double),1,fp);
    fread(&spam_gamma,sizeof(double),1,fp);
    fread(&lucy_tb_dr,sizeof(double),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&spam_flag,1,MPI_INT,0,world);
  MPI_Bcast(&flow_term_flag,1,MPI_INT,0,world);
  MPI_Bcast(&Kpd_flag,1,MPI_INT,0,world);
  MPI_Bcast(&Gpd_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&comp_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mem_stat_epsilon,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&prot_stat_epsilon,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mem_stat,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&prot_stat,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&pcov_epsilon,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&pcov_eta0,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&spam_gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&lucy_tb_dr,1,MPI_DOUBLE,0,world);

  if (comp_flag) {
    if (me == 0) {
      fread(&comp_rcut,sizeof(double),1,fp);
    }
    MPI_Bcast(&comp_rcut,1,MPI_DOUBLE,0,world);

    comp_rcutsq = comp_rcut*comp_rcut;
    comp_rcut_inv = 1.0/comp_rcut;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

/*void PairEM2::write_data(FILE *fp)
{
  int i;
  for (i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d",i);
    if (lj216_pot_flag[i][i]) fprintf(fp," %g %g %g", lj216_epsilon[i][i],lj216_k0[i][i],lj216_sigma[i][i]);
    if (lucy_pot_flag[i][i]) fprintf(fp," %g %g", lucy_epsilon[i][i],lucy_sigma[i][i]);
    if (bend_pot_flag[i][i]) fprintf(fp," %g %g %g", bend_epsilon[i][i],bend_k0[i][i],bend_gamma_epsilon[i][i]);
//    if (bend_pot_flag[i][i]) fprintf(fp," %g %g %g %g", bend_epsilon[i][i],bend_k0[i][i],bend_sigmasq[i][i],bend_gamma_epsilon[i][i]);
    if (olig_pot_flag[i][i]) fprintf(fp," %g %d", olig_epsilon[i][i],aolig[i][i]);
    if (ic_pot_flag[i]) fprintf(fp," %g %g %g", ic_lambda_m[i], ic_lambda_k[i], ic_gamma_epsilon[i]);
    if (cc_pot_flag[i]) fprintf(fp," %g %g", cc_epsilon[i], cc_zeta0[i]);
    fprintf(fp,"\n");
}*/

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

/*void PairEM2::write_data_all(FILE *fp)
{
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d",i,j);
      if (lj216_pot_flag[i][j]) fprintf(fp," %g %g %g", lj216_epsilon[i][j],lj216_k0[i][j],lj216_sigma[i][j]);
      if (lucy_pot_flag[i][j]) fprintf(fp," %g %g", lucy_epsilon[i][j],lucy_sigma[i][j]);
      if (bend_pot_flag[i][j]) fprintf(fp," %g %g %g", bend_epsilon[i][j],bend_k0[i][j],bend_gamma_epsilon[i][j]);
//      if (bend_pot_flag[i][j]) fprintf(fp," %g %g %g %g", bend_epsilon[i][j],bend_k0[i][j],bend_sigmasq[i][j],bend_gamma_epsilon[i][j]);
      if (olig_pot_flag[i][j]) fprintf(fp," %g %d", olig_epsilon[i][j],aolig[i][j]);
      if (ic_pot_flag[i]) fprintf(fp," %g %g %g", ic_lambda_m[i], ic_lambda_k[i], ic_gamma_epsilon[i]);
      if (cc_pot_flag[i]) fprintf(fp," %g %g", cc_epsilon[i], cc_zeta0[i]);
      fprintf(fp," %g\n", cut[i][j]);
}*/

int PairEM2::pack_forward_comm(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  if (comm_ind == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = rho_m[j];
      buf[m++] = rho_b[j];
      buf[m++] = n_c[j];
    }
  }
  if (comm_ind == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = xi_m[j][0];
      buf[m++] = xi_m[j][1];
      buf[m++] = xi_m[j][2];
      buf[m++] = xi_b[j][0];
      buf[m++] = xi_b[j][1];
      buf[m++] = xi_b[j][2];
    }
  }
  return m;
}

void PairEM2::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (comm_ind == 1) {
    for (i = first; i < last; i++) {
      rho_m[i] = buf[m++];
      rho_b[i] = buf[m++];
      n_c[i] = buf[m++];
    }
  }
  if (comm_ind == 2) {
    for (i = first; i < last; i++) {
      xi_m[i][0] = buf[m++];
      xi_m[i][1] = buf[m++];
      xi_m[i][2] = buf[m++];
      xi_b[i][0] = buf[m++];
      xi_b[i][1] = buf[m++];
      xi_b[i][2] = buf[m++];
    }
  }
}

int PairEM2::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (comm_ind == 1) {
    for (i = first; i < last; i++) {
      buf[m++] = rho_m[i];
      buf[m++] = rho_b[i];
      buf[m++] = n_c[i];
    }
  }
  if (comm_ind == 2) {
    for (i = first; i < last; i++) {
      buf[m++] = xi_m[i][0];
      buf[m++] = xi_m[i][1];
      buf[m++] = xi_m[i][2];
      buf[m++] = xi_b[i][0];
      buf[m++] = xi_b[i][1];
      buf[m++] = xi_b[i][2];
    }
  }
  return m;
}

void PairEM2::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (comm_ind == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      rho_m[j] += buf[m++];
      rho_b[j] += buf[m++];
      n_c[j] += buf[m++];
    }
  }
  if (comm_ind == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      xi_m[j][0] += buf[m++];
      xi_m[j][1] += buf[m++];
      xi_m[j][2] += buf[m++];
      xi_b[j][0] += buf[m++];
      xi_b[j][1] += buf[m++];
      xi_b[j][2] += buf[m++];
    }
  }
}

void *PairEM2::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"nenergy") == 0) {
    nEnergy = nEnergyTerms;
    return (void *) &nEnergy;
  }
  if (strcmp(str,"mem_stat") == 0) return (void *) &mem_stat;
  if (strcmp(str,"prot_stat") == 0) return (void *) &prot_stat;

  dim = 1;
  if (strcmp(str,"rho_m") == 0) return (void *) rho_m;
  if (strcmp(str,"rho_b") == 0) return (void *) rho_b;
  if (strcmp(str,"n_c") == 0) return (void *) n_c;
  if (strcmp(str,"energy") == 0) return (void *) energy_all;

  dim = 2;
  if (strcmp(str,"xi_m") == 0) return (void *) xi_m;
  if (strcmp(str,"xi_b") == 0) return (void *) xi_b;

  return NULL;
}
