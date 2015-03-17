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
    memory->destroy(lucy_pot_flag);
    memory->destroy(lucy_epsilon);
    memory->destroy(lucy_sigma);
    memory->destroy(lucy_sigmasq);
    memory->destroy(bend_pot_flag);
    memory->destroy(bend_epsilon);
//    memory->destroy(bend_sigmasq);
    memory->destroy(bend_k0);
    memory->destroy(bend_gamma_epsilon);
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
    memory->destroy(prot_comp_pot_flag);
    memory->destroy(prot_comp_epsilon);
    memory->destroy(prot_comp_xi_epsilon);
    memory->destroy(mem_stat_flag);
    memory->destroy(prot_stat_flag);

    memory->destroy(rho_m);
    memory->destroy(rho_b);
    memory->destroy(xi_m);
    memory->destroy(xi_b);
  }
}

/* ---------------------------------------------------------------------- */

void PairEM2::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl,one_eng,r,rsq,rinv,r2inv,r4inv,r16inv,fpair,factor_lj,epsilon,sigma,eng;
  double fforce[3],gb_force[3],torque[3],ttor[3],rtor[3],r12[3],ru12[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp;

  // Quaternions and orientations
  double *iquat,*jquat;
  double ai[3][3],aj[3][3],normi[3],normj[3],nti[3],ntj[3];

  // SPAM
  double spam_factor;
  bool b_spam_force, b_spam_force2;
  double iphi_m,jphi_m,iphi_b,jphi_b;
  double dpm,dpb,idphi[2],jdphi[2];
  double rho_one,rho_ave_inv,xi_one,dw;
  // For Stat calculations
  double sum_phi_m, sum_phi_b, sum_phi_m_all, sum_phi_b_all;
  double n_phi_m, n_phi_b, n_phi_m_all, n_phi_b_all;

  // Lucy potential
  double epsq2;

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
  double npall = newton_pair ? nall : nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // calculate per-atom densities and gradients
  // if Membrain Composition or Protein Composition potentials are on

  if (spam_flag && comp_flag) {

    // grow energy array if necessary

    if (atom->nmax > nmax) {
      memory->destroy(rho_m);
      memory->destroy(rho_b);
      memory->destroy(xi_m);
      memory->destroy(xi_b);
      nmax = atom->nmax;
      memory->create(rho_m,nmax,"pair:rho_m");
      memory->create(rho_b,nmax,"pair:rho_b");
      memory->create(xi_m,nmax,3,"pair:xi_m");
      memory->create(xi_b,nmax,3,"pair:xi_b");
    }

    // zero out densities and gradients

    for (i = 0; i < npall; i++) {
      rho_m[i] = 0.0;
      rho_b[i] = 0.0;
      xi_m[i][0] = 0.0;
      xi_m[i][1] = 0.0;
      xi_m[i][2] = 0.0;
      xi_b[i][0] = 0.0;
      xi_b[i][1] = 0.0;
      xi_b[i][2] = 0.0;
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
        factor_lj = special_lj[sbmask(j)];
        j &= NEIGHMASK;

        // r12 = center to center vector

        r12[0] = xtmp - x[j][0];
        r12[1] = ytmp - x[j][1];
        r12[2] = ztmp - x[j][2];
        rsq = MathExtra::dot3(r12,r12);
        jtype = type[j];

        // compute if less than cutoff

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
        factor_lj = special_lj[sbmask(j)];
        j &= NEIGHMASK;

        // r12 = center to center vector

        r12[0] = xtmp - x[j][0];
        r12[1] = ytmp - x[j][1];
        r12[2] = ztmp - x[j][2];
        rsq = MathExtra::dot3(r12,r12);
        jtype = type[j];

        // compute if less than cutoff

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
  }

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
        spam_factor = 0.5*(MIN(MAX(iphi_b,-1.0),0.0) + MIN(MAX(jphi_b,-1.0),0.0));
        b_spam_force = (iphi_b<=0.0 && iphi_b>=-1.0) || (newton_pair && jphi_b<=0.0 && jphi_b>=-1.0);
        b_spam_force2 = (iphi_b<=1.0 && iphi_b>=-1.0) || (newton_pair && jphi_b<=1.0 && jphi_b>=-1.0);

        // 2-16 LJ-like term for membrane-membrane interaction
        // Forces, energy, and dphi_b are calculated
        if (lj216_pot_flag[itype][jtype]) {
          epsilon = 1.0;
//          if (spam_flag && iphi_b<=0.0) epsilon *= (1 - lj216_k0[itype][jtype]*MAX(iphi_b,-1.0));
          if (spam_flag) epsilon *= (1 - lj216_k0[itype][jtype]*spam_factor);

          r2inv = 1.0/rsq;
          r4inv = r2inv*r2inv;
          r16inv = r4inv*r4inv*r4inv*r4inv;
          fpair = epsilon*(r16inv*lj1[itype][jtype]-r2inv*lj2[itype][jtype]);
          fpair *= r2inv;

          eng = (r16inv*lj3[itype][jtype]-r2inv*lj4[itype][jtype]);
          if (eflag) one_eng += epsilon*eng;
//          printf("%d %d %f\n",atom->tag[i], atom->tag[j], epsilon*eng);
          fforce[0] += r12[0]*fpair;
          fforce[1] += r12[1]*fpair;
          fforce[2] += r12[2]*fpair;

          // SPAM variable dynamics
          // Calculating -dphi/dt for phi_b only
          if (spam_flag && b_spam_force) {
            // Contribution from the derivative of prefactor
            dpb = 0.5*lj216_k0[itype][jtype]*eng;

            if (iphi_b<=0.0 && iphi_b>=-1.0) idphi[1] += dpb;
            if (jphi_b<=0.0 && jphi_b>=-1.0) jdphi[1] += dpb;
          }
        }

        // Repulsive term for solvent-solvent and membrane-solvet interaction
        // This term has form of Lucy kernel function
        // Only forces and energy are calculated
        if (lucy_pot_flag[itype][jtype] && r<lucy_sigma[itype][jtype]) {
          epsilon = lucy_epsilon[itype][jtype];

          sigma = r/lucy_sigma[itype][jtype];
          epsq2 = epsilon*(1-sigma)*(1-sigma);
          fpair = 12*epsq2/lucy_sigmasq[itype][jtype];

          if (eflag) one_eng += epsq2*(1-sigma)*(1+3*sigma);
//          printf("%d %d %f\n",atom->tag[i], atom->tag[j], epsq2*(1-sigma)*(1+3*sigma));
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
          if (spam_flag) epsilon_sigmasq_r2inv *= (1 - bend_k0[itype][jtype]*spam_factor);

//          if (spam_flag) epsilon_sigmasq_r2inv *= (1 - 0.5*bend_k0[itype][jtype]*(MIN(MAX(iphi_b,-1.0),0.0) + MIN(MAX(jphi_b,-1.0),0.0)));

//          epsilon = 4.0*bend_epsilon[itype][jtype];
//          if (spam_flag) epsilon *= (1 - 0.5*bend_k0[itype][jtype]*(MIN(MAX(iphi_b,-1.0),0.0) + MIN(MAX(jphi_b,-1.0),0.0)));
//          epsilon_sigmasq_r2inv = epsilon*(bend_sigmasq[itype][jtype]*r2inv);

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
            fi = MAX(MIN(-iphi_b,1.0),0.0);
            fj = MAX(MIN(-jphi_b,1.0),0.0);
            gamma_epsilon *= 0.25*(fi + fj);
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

//          fpair = -uij*r2inv;
//          fpair -= (thetai - gamma_r)*(thetai*r2inv + gamma*rinv - 2*gamma_epsilon*alpha_sq*rinv);
//          fpair -= (thetaj + gamma_r)*(thetaj*r2inv - gamma*rinv + 2*gamma_epsilon*alpha_sq*rinv);

          if (eflag) one_eng += epsilon_sigmasq_r2inv*uij;
//          printf("%d %d %f\n",atom->tag[i], atom->tag[j], epsilon_sigmasq_r2inv*uij);

          // Adding all forces
          fforce[0] += fpair*r12[0] + gb_force[0];
          fforce[1] += fpair*r12[1] + gb_force[1];
          fforce[2] += fpair*r12[2] + gb_force[2];
//          fforce[0] += fpair*r12[0] + fthi*thetai[0] + fthj*thetaj[0] + fali*alphai[0] + falj*alphaj[0];
//          fforce[1] += fpair*r12[1] + fthi*thetai[1] + fthj*thetaj[1] + fali*alphai[1] + falj*alphaj[1];
//          fforce[2] += fpair*r12[2] + fthi*thetai[2] + fthj*thetaj[2] + fali*alphai[2] + falj*alphaj[2];

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

            // Contribution from the derivative of prefactor
            dpb = 0.5*bend_epsilon[itype][jtype]*bend_k0[itype][jtype]*r2inv*uij;
            // Contribution from the derivative of gamma
            dpb -= 0.25*epsilon_sigmasq_r2inv*bend_gamma_epsilon[itype][jtype]*(thetai_gamma_r - thetaj_gamma_r)*alpha_sq*r;
            
            if (iphi_b<=0.0 && iphi_b>=-1.0) idphi[1] += dpb;
            if (jphi_b<=0.0 && jphi_b>=-1.0) jdphi[1] += dpb;
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
            
            if (eflag) one_eng -= epsilon*MathExtra::dot3(nti, ntj)/((double)aolig[itype][jtype]);

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

          epsilon = ic_lambda_k[itype]*MIN(MAX(iphi_b,-1.0),1.0) + ic_lambda_k[jtype]*MIN(MAX(jphi_b,-1.0),1.0);

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
          fi = MAX(MIN(-iphi_b,1.0),0.0);
          fj = MAX(MIN(-jphi_b,1.0),0.0);
          gamma_epsilon = 0.125*(ic_gamma_epsilon[itype]*fi + ic_gamma_epsilon[jtype]*fj);

          gamma = gamma_epsilon*alpha_sq;
          gamma_r = gamma*r;

          thetai_gamma_r = thetai - gamma_r;
          thetaj_gamma_r = thetaj + gamma_r;

          uij = thetai_gamma_r*thetai_gamma_r + thetaj_gamma_r*thetaj_gamma_r;

          // Spherical force
          fpair = -(thetai_gamma_r*thetai + thetaj_gamma_r*thetaj)*r2inv;
          fpair -= (thetai_gamma_r - thetaj_gamma_r)*(gamma - 2.0*gamma_epsilon*alpha_sq)*rinv;
          fpair *= 2.0*epsilon;

          fthi = thetai_gamma_r*rinv;
          fthj = thetaj_gamma_r*rinv;

          fali = -2.0*gamma_epsilon*(thetai_gamma_r - thetaj_gamma_r)*alphai;
          falj = -2.0*gamma_epsilon*(thetai_gamma_r - thetaj_gamma_r)*alphaj;

          // Aspherical forces
          gb_force[0] = fthi*normi[0] + fthj*normj[0] + fali*nti[0] + falj*ntj[0];
          gb_force[1] = fthi*normi[1] + fthj*normj[1] + fali*nti[1] + falj*ntj[1];
          gb_force[2] = fthi*normi[2] + fthj*normj[2] + fali*nti[2] + falj*ntj[2];

          gb_force[0] *= 2.0*epsilon;
          gb_force[1] *= 2.0*epsilon;
          gb_force[2] *= 2.0*epsilon;

          if (eflag) one_eng -= epsilon*uij;
//          printf("%d %d %f\n",atom->tag[i], atom->tag[j], -epsilon*uij);

          // Adding all forces
          fforce[0] += fpair*r12[0] + gb_force[0];
          fforce[1] += fpair*r12[1] + gb_force[1];
          fforce[2] += fpair*r12[2] + gb_force[2];

          // Calculate torque on particle i

          // Torque from derivative over (ru, normi)
          tor_thi = -2.0*epsilon*thetai_gamma_r;
          MathExtra::cross3(ru12, normi, torque);
          ttor[0] += tor_thi*torque[0];
          ttor[1] += tor_thi*torque[1];
          ttor[2] += tor_thi*torque[2];

          // Torque from derivative over (ru, nti)
          tor_ali = -2.0*epsilon*fali*r;
          MathExtra::cross3(ru12, nti, torque);
          ttor[0] += tor_ali*torque[0];
          ttor[1] += tor_ali*torque[1];
          ttor[2] += tor_ali*torque[2];

          // Calculate torque on particle j

          // Torque from derivative over (ru, normj)
          tor_thj = -2.0*epsilon*thetaj_gamma_r;
          MathExtra::cross3(ru12, normj, torque);
          rtor[0] += tor_thj*torque[0];
          rtor[1] += tor_thj*torque[1];
          rtor[2] += tor_thj*torque[2];

          // Torque from derivative over (ru, ntj)
          tor_alj = -2.0*epsilon*falj*r;
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
            dpb = 0.25*epsilon*(thetai_gamma_r - thetaj_gamma_r)*alpha_sq*r;
            if (iphi_b<=1.0 && iphi_b>=-1.0) {
              idphi[1] += ic_lambda_k[itype]*uij;
              if (iphi_b<=0.0) idphi[1] += ic_gamma_epsilon[itype]*dpb;
            }
            if (newton_pair || j < nlocal) {
              if (jphi_b<=1.0 && jphi_b>=-1.0) {
                jdphi[1] += ic_lambda_k[jtype]*uij;
                if (jphi_b<=0.0) jdphi[1] += ic_gamma_epsilon[jtype]*dpb;
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

/*          if (atom->tag[i]==1) {
            printf("tagi=%d tagj=%d irho=%f jrho=%f ixi={%f %f %f} jxi={%f %f %f} r=%f r12={%f %f %f} sigma=%f dw=%f xi_epsilon=%f dphi=%f\n", atom->tag[i], atom->tag[j], rho_m[i], rho_m[j], xi_m[i][0], xi_m[i][1], xi_m[i][2], xi_m[j][0], xi_m[j][1], xi_m[j][2], r, r12[0], r12[1], r12[2], sigma, dw, mem_comp_xi_epsilon[itype]*mass[jtype], mem_comp_xi_epsilon[itype]*mass[jtype]*dpm);
          }
          if (atom->tag[j]==1) {
            printf("tagi=%d tagj=%d irho=%f jrho=%f ixi={%f %f %f} jxi={%f %f %f} r=%f r12={%f %f %f} sigma=%f dw=%f xi_epsilon=%f dphi=%f\n", atom->tag[i], atom->tag[j], rho_m[i], rho_m[j], xi_m[i][0], xi_m[i][1], xi_m[i][2], xi_m[j][0], xi_m[j][1], xi_m[j][2], r, r12[0], r12[1], r12[2], sigma, dw, mem_comp_xi_epsilon[jtype]*mass[itype], mem_comp_xi_epsilon[jtype]*mass[itype]*dpm);
          }*/
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

/*          if (atom->tag[i]==105) {
            printf("tagi=%d tagj=%d irho=%f jrho=%f ixi={%f %f %f} jxi={%f %f %f} r=%f r12={%f %f %f} sigma=%f dw=%.12f xi_epsilon=%f dphi=%.12f\n", atom->tag[i], atom->tag[j], rho_b[i], rho_b[j], xi_b[i][0], xi_b[i][1], xi_b[i][2], xi_b[j][0], xi_b[j][1], xi_b[j][2], r, r12[0], r12[1], r12[2], sigma, dw, prot_comp_xi_epsilon[itype]*mass[jtype], prot_comp_xi_epsilon[itype]*mass[jtype]*dpb);
          }
          if (atom->tag[j]==105) {
            printf("tagi=%d tagj=%d irho=%f jrho=%f ixi={%f %f %f} jxi={%f %f %f} r=%f r12={%f %f %f} sigma=%f dw=%.12f xi_epsilon=%f dphi=%.12f\n", atom->tag[i], atom->tag[j], rho_b[i], rho_b[j], xi_b[i][0], xi_b[i][1], xi_b[i][2], xi_b[j][0], xi_b[j][1], xi_b[j][2], r, r12[0], r12[1], r12[2], sigma, dw, prot_comp_xi_epsilon[jtype]*mass[itype], prot_comp_xi_epsilon[jtype]*mass[itype]*dpb);
          }*/
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
      if (eflag) one_eng += ic_lambda_m[itype]*MIN(MAX(iphi_b,-1.0),1.0);
//      printf("itag=%d iphi_b=%f/%f eng=%f\n", atom->tag[i], iphi_b, MIN(MAX(iphi_b,-1.0),1.0), ic_lambda_m[itype]*MIN(MAX(iphi_b,-1.0),1.0));

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

      if (eflag) one_eng -= cc_epsilon[itype]*zeta_m*zeta_b;

      idphi[0] += dpm;
      idphi[1] += dpb;
    }

    // Membrain composition potential
    // phi_m^10/10
    // dphi_m from one-body well potential and total energy are calculated
    // dphi_m from gradient is calculated in the double loop above
    if (spam_flag && mem_comp_pot_flag[itype]) {
      dpm = mem_comp_epsilon[itype]*pow(iphi_m,9.0);

      // phi_m gradient term energy
      if (eflag) one_eng += mem_comp_xi_epsilon[itype]*MathExtra::dot3(xi_m[i],xi_m[i]);
      // Well potential
      if (eflag) one_eng += mem_comp_epsilon[itype]*0.1*pow(iphi_m,10.0);

//      if (atom->tag[i]==2727) printf("%d %d %f %f %f\n",update->ntimestep,atom->tag[i],iphi_m,dpm,mem_comp_epsilon[itype]*0.1*pow(iphi_m,10.0));

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
//      dpb = prot_comp_epsilon[itype]*(6.0*pow(iphi_b,5.0) + 2.0*iphi_b);

      // phi_b gradient term energy
      if (eflag) one_eng += prot_comp_xi_epsilon[itype]*MathExtra::dot3(xi_b[i],xi_b[i]);
      // Well potential
      if (eflag) one_eng += prot_comp_epsilon[itype]*phi2*(phi4 + 1.0);
//      if (eflag) one_eng += prot_comp_epsilon[itype]*(pow(iphi_b,6.0) + pow(iphi_b,2.0));

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

  // Propogate mem_stat and prot_stat
  // stat += epsilon*sum_phi/n_phi

  if (n_phi_m_all>0.0) mem_stat += update->dt*mem_stat_epsilon*sum_phi_m_all/n_phi_m_all;
  if (n_phi_b_all>0.0) prot_stat += update->dt*prot_stat_epsilon*sum_phi_b_all/n_phi_b_all;

//  if(comm->me==0) printf("MEM  STAT step: %d sum_phi: %f n_phi: %f dstat: %f stat: %f\n", update->ntimestep, sum_phi_m_all, n_phi_m_all, update->dt*mem_stat_epsilon*sum_phi_m_all/n_phi_m_all, mem_stat);
//  if(comm->me==0) printf("PROT STAT step: %d sum_phi: %f n_phi: %f dstat: %f stat: %f\n", update->ntimestep, sum_phi_b_all, n_phi_b_all, update->dt*prot_stat_epsilon*sum_phi_b_all/n_phi_b_all, prot_stat);

  // Second one-body loop
  // idphi here is NOT multipied by spam_gamma
  for (i=0; i < nlocal; i++) {
    itype = type[i];

    idphi[0] = idphi[1] = 0.0;

    // Subtract stat from dphi
    if (spam_flag) {
      
      idphi[0] -= mem_stat;
      idphi[1] -= prot_stat;
    }

    // Velocity times composition gradient
    // Possibly need to move to integration
    if (spam_flag && flow_term_flag) {
//      if (atom->tag[i]==2727) printf("FLOW %d %d %f\n", update->ntimestep, atom->tag[i], MathExtra::dot3(v[i],xi_m[i]));
      if (mem_comp_pot_flag[itype]) idphi[0] += MathExtra::dot3(v[i],xi_m[i]);
      if (prot_comp_pot_flag[itype]) idphi[1] += MathExtra::dot3(v[i],xi_b[i]);
    }

    dphi[i][0] += idphi[0];
    dphi[i][1] += idphi[1];
  }

/*  for (i=0;i<nall;i++) {
    if (atom->tag[i]==105) {
      printf("%d %f %f\n", atom->tag[i], dphi[i][0], dphi[i][1]);
    }
  }*/

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
  memory->create(lucy_pot_flag,n+1,n+1,"pair:lucy_pot_flag");
  memory->create(lucy_epsilon,n+1,n+1,"pair:lucy_epsilon");
  memory->create(lucy_sigma,n+1,n+1,"pair:lucy_sigma");
  memory->create(lucy_sigmasq,n+1,n+1,"pair:lucy_sigmasq");
  memory->create(bend_pot_flag,n+1,n+1,"pair:bend_pot_flag");
  memory->create(bend_epsilon,n+1,n+1,"pair:bend_epsilon");
//  memory->create(bend_sigmasq,n+1,n+1,"pair:bend_sigmasq");
  memory->create(bend_k0,n+1,n+1,"pair:bend_k0");
  memory->create(bend_gamma_epsilon,n+1,n+1,"pair:");
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
  memory->create(prot_comp_pot_flag,n+1,"pair:prot_comp_pot_flag");
  memory->create(prot_comp_epsilon,n+1,"pair:prot_comp_epsilon");
  memory->create(prot_comp_xi_epsilon,n+1,"pair:prot_comp_xi_epsilon");
  memory->create(mem_stat_flag,n+1,"pair:mem_stat_flag");
  memory->create(prot_stat_flag,n+1,"pair:prot_stat_flag");

  comp_flag = 0;
  flow_term_flag = 0;
  for (int i = 1; i <= n; i++) {
    ic_pot_flag[i] = 0;
    cc_pot_flag[i] = 0;
    mem_comp_pot_flag[i] = 0;
    prot_comp_pot_flag[i] = 0;
    mem_stat_flag[i] = 0;
    prot_stat_flag[i] = 0;
    for (int j = 1; j <= n; j++) {
      lj216_pot_flag[i][j] = 0;
      lucy_pot_flag[i][j] = 0;
      bend_pot_flag[i][j] = 0;
      olig_pot_flag[i][j] = 0;
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

  // Need to check if this is a right place for this
  mem_stat = 0.0;
  prot_stat = 0.0;

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
}

/* ----------------------------------------------------------------------*/

void PairEM2::read_parameters()
{
  read_par_flag = false;

  int i, j;
  int ilo,ihi,jlo,jhi;
  double epsilon_val, epsilon2_val, sigma_val, k0_val, geps_val;
  double lambda_m_val, lambda_k_val, zeta0_val;
  int aolig_val;
  int flag1, flag2;

  FILE *file;
  int file_state, narg=0, iarg_line=-1;
  char ln[1024], *line, *arg[16];
  enum File_States{FS_NONE=0, FS_LJ216, FS_LUCY, FS_BEND, FS_OLIG, FS_COUP_IC, FS_COUP_CC, FS_MEM_COMP, FS_BAR_COMP, FS_COMP_CUT, FS_COMP_STAT, FS_FLOW, FS_SPAM};

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
//    if (narg<2) error->all(FLERR,"Pair_style EM2: Error reading parameter file!");

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
//          bend_sigmasq[i][j] = sigma_val*sigma_val;
          bend_gamma_epsilon[i][j] = geps_val;
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
      if (narg!=3) error->all(FLERR,"Pair_style EM2: Wrong format in coefficient file (Membrain Composition Pot)");
      if (me==0) print_log("Membrain Composition potential flag on\n");
      epsilon_val = atof(arg[1]);
      epsilon2_val = atof(arg[2]);

      force->bounds(arg[0],atom->ntypes,ilo,ihi);

      for (i = ilo; i <= ihi; i++) {
        mem_comp_pot_flag[i] = 1;
        mem_comp_xi_epsilon[i] = epsilon_val;
        mem_comp_epsilon[i] = epsilon2_val;
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
    case FS_NONE:
      iarg_line = -1;
      if (strcmp(line, "[Lennard-Jones_2-16]")==0)
        file_state = FS_LJ216;
      else if (strcmp(line, "[Lucy_Excluded_Volume]")==0)
        file_state = FS_LUCY;
      else if (strcmp(line, "[EM2_Bending]")==0)
        file_state = FS_BEND;
      else if (strcmp(line, "[Oligomerization_Energy]")==0)
        file_state = FS_OLIG;
      else if (strcmp(line, "[Intrinsic_Curvature_Coupling]")==0)
        file_state = FS_COUP_IC;
      else if (strcmp(line, "[Composition_Coupling]")==0)
        file_state = FS_COUP_CC;
      else if (strcmp(line, "[Membrane_Composition]")==0)
        file_state = FS_MEM_COMP;
      else if (strcmp(line, "[Protein_Composition]")==0)
        file_state = FS_BAR_COMP;
      else if (strcmp(line, "[Composition_Cutoff]")==0)
        file_state = FS_COMP_CUT;
      else if (strcmp(line, "[Composition_Stat]")==0)
        file_state = FS_COMP_STAT;
      else if (strcmp(line, "[Flow_Term_Flag]")==0)
        file_state = FS_FLOW;
      else if (strcmp(line, "[SPAM_Flag]")==0)
        file_state = FS_SPAM;
      break;
    }
  }
  fclose(file);
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

  if (flag && !comp_flag) error->all(FLERR,"Pair em2 requires composition cutoff if any corresponding potentials are on");
  if (!flag) comp_flag = 0;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEM2::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

//  if (setflag[i][j] == 0) cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  if (lj216_pot_flag[i][j]) {
/*    if (setflag[i][j] == 0) {
      lj216_pot_flag[i][j] = 1;
      lj216_epsilon[i][j] = mix_energy(lj216_epsilon[i][i],lj216_epsilon[j][j],
                                 lj216_sigma[i][i],lj216_sigma[j][j]);
      lj216_sigma[i][j] = mix_distance(lj216_sigma[i][i],lj216_sigma[j][j]);
    }*/

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

  if (lucy_pot_flag[i][j]) {
/*    if (setflag[i][j] == 0) {
      lucy_epsilon[i][j] = mix_energy(lucy_epsilon[i][i],lucy_epsilon[j][j],
                                 lucy_sigma[i][i],lucy_sigma[j][j]);
      lucy_sigma[i][j] = mix_distance(lucy_sigma[i][i],lucy_sigma[j][j]);
    }*/

    if (lucy_sigma[i][j]>cut[i][j]) error->all(FLERR,"Lucy Excluded Volume cutoff cannot be larger than overal cutoff for a pair type");

    lucy_sigmasq[i][j] = pow(lucy_sigma[i][j], 2.0);

    lucy_pot_flag[j][i] = lucy_pot_flag[i][j];
    lucy_epsilon[j][i] = lucy_epsilon[i][j];
    lucy_sigma[j][i] = lucy_sigma[i][j];
    lucy_sigmasq[j][i] = lucy_sigmasq[i][j];
  }

  if (bend_pot_flag[i][j]) {
/*    if (setflag[i][j] == 0) {
      bend_epsilon[i][j] = mix_energy(bend_epsilon[i][i],bend_epsilon[j][j],
                                 bend_sigma[i][i],bend_sigma[j][j]);
      bend_sigma[i][j] = mix_distance(bend_sigma[i][i],bend_sigma[j][j]);
      bend_k0[i][j] = mix_energy(bend_k0[i][i], bend_k0[j][j], 1.0, 1.0);
      bend_gamma_epsilon[i][j] = mix_energy(bend_gamma_epsilon[i][i], bend_gamma_epsilon[j][j], 1.0, 1.0);
    }*/
   
    bend_pot_flag[j][i] = bend_pot_flag[i][j]; 
    bend_epsilon[j][i] = bend_epsilon[i][j];
//    bend_sigmasq[j][i] = bend_sigmasq[i][j];
    bend_k0[j][i] = bend_k0[i][j];
    bend_gamma_epsilon[j][i] = bend_gamma_epsilon[i][j];
  }

  if (olig_pot_flag[i][j]) {
/*    if (setflag[i][j] == 0) {
      olig_epsilon[i][j] = mix_energy(olig_epsilon[i][i], olig_epsilon[j][j], 1.0, 1.0);
      aolig[i][j] = 1;
    }*/

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
    fwrite(&prot_comp_pot_flag[i],sizeof(int),1,fp);
    fwrite(&mem_stat_flag[i],sizeof(int),1,fp);
    fwrite(&prot_stat_flag[i],sizeof(int),1,fp);
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
      fwrite(&lucy_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&bend_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&olig_pot_flag[i][j],sizeof(int),1,fp);
      fwrite(&cut[i][j],sizeof(double),1,fp);
      if (lj216_pot_flag[i][j]) {
        fwrite(&lj216_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&lj216_k0[i][j],sizeof(double),1,fp);
        fwrite(&lj216_sigma[i][j],sizeof(double),1,fp);
      }
      if (lucy_pot_flag[i][j]) {
        fwrite(&lucy_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&lucy_sigma[i][j],sizeof(double),1,fp);
      }
      if (bend_pot_flag[i][j]) {
        fwrite(&bend_epsilon[i][j],sizeof(double),1,fp);
        fwrite(&bend_k0[i][j],sizeof(double),1,fp);
//        fwrite(&bend_sigmasq[i][j],sizeof(double),1,fp);
        fwrite(&bend_gamma_epsilon[i][j],sizeof(double),1,fp);
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
  allocate();

  int i,j;
  int me = comm->me;

  // Single loop
  for (i = 1; i <= atom->ntypes; i++) {
    if (me == 0) {
      fread(&ic_pot_flag[i],sizeof(int),1,fp);
      fread(&cc_pot_flag[i],sizeof(int),1,fp);
      fread(&mem_comp_pot_flag[i],sizeof(int),1,fp);
      fread(&prot_comp_pot_flag[i],sizeof(int),1,fp);
      fread(&mem_stat_flag[i],sizeof(int),1,fp);
      fread(&prot_stat_flag[i],sizeof(int),1,fp);
    }
    MPI_Bcast(&ic_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&cc_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&mem_comp_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&prot_comp_pot_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&mem_stat_flag[i],1,MPI_INT,0,world);
    MPI_Bcast(&prot_stat_flag[i],1,MPI_INT,0,world);

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
      }
      MPI_Bcast(&mem_comp_epsilon[i],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&mem_comp_xi_epsilon[i],1,MPI_DOUBLE,0,world);
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
        fread(&lucy_pot_flag[i][j],sizeof(int),1,fp);
        fread(&bend_pot_flag[i][j],sizeof(int),1,fp);
        fread(&olig_pot_flag[i][j],sizeof(int),1,fp);
      }
      MPI_Bcast(&lj216_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&lucy_pot_flag[i][j],1,MPI_INT,0,world);
      MPI_Bcast(&bend_pot_flag[i][j],1,MPI_INT,0,world);
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

      if (lucy_pot_flag[i][j]) {
        if (me == 0) {
          fread(&lucy_epsilon[i][j],sizeof(double),1,fp);
          fread(&lucy_sigma[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&lucy_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lucy_sigma[i][j],1,MPI_DOUBLE,0,world);
      }

      if (bend_pot_flag[i][j]) {
        if (me == 0) {
          fread(&bend_epsilon[i][j],sizeof(double),1,fp);
          fread(&bend_k0[i][j],sizeof(double),1,fp);
//          fread(&bend_sigmasq[i][j],sizeof(double),1,fp);
          fread(&bend_gamma_epsilon[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&bend_epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&bend_k0[i][j],1,MPI_DOUBLE,0,world);
//        MPI_Bcast(&bend_sigmasq[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&bend_gamma_epsilon[i][j],1,MPI_DOUBLE,0,world);
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
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&comp_flag,sizeof(int),1,fp);
  fwrite(&mem_stat_epsilon,sizeof(double),1,fp);
  fwrite(&prot_stat_epsilon,sizeof(double),1,fp);
  fwrite(&mem_stat,sizeof(double),1,fp);
  fwrite(&prot_stat,sizeof(double),1,fp);
  fwrite(&spam_gamma,sizeof(double),1,fp);

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
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&comp_flag,sizeof(int),1,fp);
    fread(&mem_stat_epsilon,sizeof(double),1,fp);
    fread(&prot_stat_epsilon,sizeof(double),1,fp);
    fread(&mem_stat,sizeof(double),1,fp);
    fread(&prot_stat,sizeof(double),1,fp);
    fread(&spam_gamma,sizeof(double),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&spam_flag,1,MPI_INT,0,world);
  MPI_Bcast(&flow_term_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&comp_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mem_stat_epsilon,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&prot_stat_epsilon,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mem_stat,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&prot_stat,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&spam_gamma,1,MPI_DOUBLE,0,world);

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
  dim = 1;
  if (strcmp(str,"rho_m") == 0) return (void *) rho_m;
  if (strcmp(str,"rho_b") == 0) return (void *) rho_b;

  dim = 2;
  if (strcmp(str,"xi_m") == 0) return (void *) xi_m;
  if (strcmp(str,"xi_b") == 0) return (void *) xi_b;

  return NULL;
}
