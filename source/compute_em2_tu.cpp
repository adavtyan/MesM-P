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

#include "compute_em2_tu.h"
#include "atom_vec_em2.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "comm.h"
#include "modify.h"
#include "pair.h"
#include "force.h"
#include "string.h"
#include "stdlib.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEM2TU::ComputeEM2TU(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4 && narg != 5) error->all(FLERR,"Illegal compute em2_tu command");

  scalar_flag = 1;
  extscalar = 0;
  
  double rcut = atof(arg[3]);
  if (rcut < 0.0) error->all(FLERR,"Illegal compute em2_tu command");
  rcutsq = rcut*rcut;
  
  plist = PL_HALF;
  if (narg>=5) {
    if (strcmp(arg[4],"comm")==0) plist = PL_HALF;
    else if (strcmp(arg[4],"no_comm")==0) plist = PL_FULL;
    else error->all(FLERR,"Illegal argument for compute em2_tu command");
  }

  nmax = 0;
  tu = NULL;
  nc = NULL;
  
  norm0[0] = 0.0;
  norm0[1] = 0.0;
  norm0[2] = 1.0;

  if (plist==PL_HALF) {
    comm_forward = 2;
    comm_reverse = 2;
  }
}

/* ---------------------------------------------------------------------- */

ComputeEM2TU::~ComputeEM2TU()
{
  memory->destroy(tu);
  memory->destroy(nc);
}

/* ---------------------------------------------------------------------- */

void ComputeEM2TU::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"em2_tu") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute em2_tu");
  
  if (force->pair == NULL)
    error->all(FLERR,"Compute em2_tu requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute em2_tu cutoff is longer than pairwise cutoff");

  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) error->all(FLERR,"Compute em2_tu requires atom style em2");

  // need an occasional half neighbor list
  
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  if (plist!=PL_HALF) {
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeEM2TU::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeEM2TU::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  
  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list);
  
  // grow per-atom arrays if needed
  
  if (atom->nmax > nmax) {
    memory->destroy(tu);
    memory->destroy(nc);
    nmax = atom->nmax;
    memory->create(tu,nmax,"compute:tu");
    memory->create(nc,nmax,"compute:nc");
  }
  
  // Either use half pair list and communication
  // or use full pair list to compute nc and tu
  
  if (plist==PL_HALF) return compute_whalf_pl();
  return compute_wfull_pl();
}

double ComputeEM2TU::compute_whalf_pl()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double tu_one, tu_sum, quat_dot;
  double *iquat, *jquat;
  double ai[3][3],aj[3][3],normi[3],normj[3];

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  
  int npall = newton_pair ? nall : nlocal;
  
  for (i=0; i < npall; i++) {
    if ((mask[i] & groupbit))
      tu[i] = 0.0;
      nc[i] = 0.0;
    }
  }
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    iquat = bonus[i].quat;
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    // Normal vector for particle i
    MathExtra::quat_to_mat_trans(iquat, ai);
    MathExtra::vecmat(norm0,ai,normi);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq<rcutsq) {
        jquat = bonus[j].quat;

        // Normal vector for particle i
        MathExtra::quat_to_mat_trans(jquat, aj);
        MathExtra::vecmat(norm0,aj,normj);
        
        quat_dot = MathExtra::dot3(normi, normj);
        tu_one = 1.5*quat_dot*quat_dot - 0.5;
        
        tu[i] += tu_one;
        nc[i] += 1.0;
        if (newton_pair || j < nlocal) {
          tu[j] += tu_one;
          nc[j] += 1.0;
        }
      }
    }
  }
  
  if (newton_pair) comm->reverse_comm_compute(this);
  comm->forward_comm_compute(this);
  
  tu_one = 0.0;
  for (i=0; i < nlocal; i++) {
    if ((mask[i] & groupbit)) {
      if (nc[i]>0.0) tu[i] /= nc[i];
      else tu[i] = 0.0;
      tu_one += tu[i];
    }
  }
  
  MPI_Allreduce(&tu_one,&tu_sum,1,MPI_DOUBLE,MPI_SUM,world);
  
  return tu_sum;
}

double ComputeEM2TU::compute_wfull_pl()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double tu_one, tu_sum, quat_dot;
  double *iquat, *jquat;
  double ai[3][3],aj[3][3],normi[3],normj[3];

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (i=0; i < nlocal; i++) {
    if ((mask[i] & groupbit))
      tu[i] = 0.0;
      nc[i] = 0.0;
    }
  }
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    iquat = bonus[i].quat;
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    // Normal vector for particle i
    MathExtra::quat_to_mat_trans(iquat, ai);
    MathExtra::vecmat(norm0,ai,normi);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq<rcutsq) {
        jquat = bonus[j].quat;

        // Normal vector for particle i
        MathExtra::quat_to_mat_trans(jquat, aj);
        MathExtra::vecmat(norm0,aj,normj);
        
        quat_dot = MathExtra::dot3(normi, normj);
        tu_one = 1.5*quat_dot*quat_dot - 0.5;
        
        tu[i] += tu_one;
        nc[i] += 1.0;
      }
    }
  }
   
  tu_one = 0.0;
  for (i=0; i < nlocal; i++) {
    if ((mask[i] & groupbit)) {
      if (nc[i]>0.0) tu[i] /= nc[i];
      else tu[i] = 0.0;
      tu_one += tu[i];
    }
  }
  
  MPI_Allreduce(&tu_one,&tu_sum,1,MPI_DOUBLE,MPI_SUM,world);
  
  return tu_sum;
}

int ComputeEM2TU::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = nc[j];
    buf[m++] = tu[j];
  }
   return m;
}

void ComputeEM2TU::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    nc[i] = buf[m++];
    tu[i] = buf[m++];
  }
}

int ComputeEM2TU::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = nc[i];
    buf[m++] = tu[i];
  }
  return m;
}

void ComputeEM2TU::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    nc[j] += buf[m++];
    tu[j] += buf[m++];
  }
}

double ComputeEM2TU::memory_usage()
{
  double bytes = 2*nmax * sizeof(double);
  return bytes;
}

