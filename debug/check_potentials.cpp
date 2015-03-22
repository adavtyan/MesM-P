#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double L = 3600.0;
double L_half = 1800.0;
double phi_g = 4.0;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

double dot(double *a, double *b);
double dist_sq(double *x1, double *x2);
double delta(double *x1, double *x2, double *dx);

class Atoms 
{
public:
  Atoms(int m, int n);
  ~Atoms();
  int n_atom;
  int n_mem;
  int *type;
  double **x;
  double **quat;
  double **phi;
};

class PairList
{
public:
  PairList();
  ~PairList();
  int make_pair_list(double rcutsq, Atoms &atoms, int n);
  int npair;
  int natom;
  int *list;
};

static PairList PVAL;

Atoms::Atoms(int m, int n)
{
  type = (int*) malloc (n * sizeof(int));
  quat = (double **)malloc(n * sizeof(double *));
  phi = (double **)malloc(n * sizeof(double *));
  x = (double **)malloc(n * sizeof(double *));
  for (int i = 0; i < n; i++) {
    quat[i] = (double *)malloc(4 * sizeof(double));
    phi[i] = (double *)malloc(2 * sizeof(double));
    x[i] = (double *)malloc(3 * sizeof(double));
  }

  n_atom = n;
  n_mem = m;
}

Atoms::~Atoms()
{
  for (int i = 0; i < n_atom; i++) {
    free(quat[i]);
    free(phi[i]);
    free(x[i]);
  }
  free(type);
  free(quat);
  free(phi);
  free(x);
}

PairList::PairList()
{
  npair = 0;
  list = NULL;
}

PairList::~PairList()
{
  if (list) free(list);
}

int PairList::make_pair_list(double rcutsq, Atoms &a, int n)
{
  int i,j;
  int ncur = 0;
  int nblock = 1000;
  double rsq;

  npair = 0;
  natom = n;
  
  ncur = nblock;
  list = (int *)malloc(2*ncur * sizeof(int));

  double **x = a.x;

  for (i=0;i<n;++i) {
    for (j=i+1;j<n;j++) {
      rsq = dist_sq(x[i], x[j]);
      if (rsq<rcutsq) {
        if (npair==ncur) {
          ncur += nblock;
          list = (int *)realloc(list, 2*ncur * sizeof(int));
        }

        list[2*npair] = i;
        list[2*npair+1] = j;
        npair++;
      }
    }
  }

  if (npair!=ncur) list = (int *)realloc(list, 2*npair * sizeof(int));

  return npair;
}

double dot(double *a, double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double dist_sq(double *x1, double *x2)
{
  double dx[3];
  dx[0] = x1[0] - x2[0];
  dx[1] = x1[1] - x2[1];
  dx[2] = x1[2] - x2[2];

  if (dx[0]>L_half) dx[0] -= L;
  else if (dx[0]<-L_half) dx[0] += L;
  if (dx[1]>L_half) dx[1] -= L;
  else if (dx[1]<-L_half) dx[1] += L;
  if (dx[2]>L_half) dx[2] -= L;
  else if (dx[2]<-L_half) dx[2] += L;

  return dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
}

double delta(double *x1, double *x2, double *dx)
{
  dx[0] = x1[0] - x2[0];
  dx[1] = x1[1] - x2[1];
  dx[2] = x1[2] - x2[2];

  if (dx[0]>L_half) dx[0] -= L;
  else if (dx[0]<-L_half) dx[0] += L;
  if (dx[1]>L_half) dx[1] -= L;
  else if (dx[1]<-L_half) dx[1] += L;
  if (dx[2]>L_half) dx[2] -= L;
  else if (dx[2]<-L_half) dx[2] += L;
  
  return dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
}

double lucy_pot(double *par, Atoms &atoms, PairList &pair)
{
  int ii,i,j;
  double r,rsq,sigma,rsq_cut,rcut_inv;
  double eng = 0.0, eng_one;

  rsq_cut = par[2]*par[2];
  rcut_inv = 1.0/par[2];

  int n = pair.npair;
  int *list = pair.list;
  int *type = atoms.type;
  double **x = atoms.x;

  for (ii=0;ii<n;++ii) {
    i = list[2*ii];
    j = list[2*ii+1];
    if (type[i]==1 && type[j]==1) continue;
    rsq = dist_sq(x[i], x[j]);
    if (rsq<rsq_cut) {
      r = sqrt(rsq);
      sigma = r*rcut_inv;
      if (type[i]==2 && type[j]==2)
        eng_one = par[1]*(1-sigma)*(1-sigma)*(1-sigma)*(1+3*sigma);
      else
        eng_one = par[0]*(1-sigma)*(1-sigma)*(1-sigma)*(1+3*sigma);
//      printf("%d %d %f\n", i+1, j+1, eng_one);
      eng += eng_one;
    }
  }

  return eng;
}

double lj216_pot(double *par, Atoms &atoms, PairList &pair)
{
  int ii,i,j;
  double r,rsq,sigma,rsq_cut,f12;
  double eng_one, eng = 0.0;

  rsq_cut = par[3]*par[3];

  int n = pair.npair;
  int *list = pair.list;
  int *type = atoms.type;
  double **x = atoms.x;
  double **phi = atoms.phi;

  for (ii=0;ii<n;++ii) {
    i = list[2*ii];
    j = list[2*ii+1];
    if (type[i]==1 && type[j]==1) {
      rsq = dist_sq(x[i], x[j]);
      if (rsq<rsq_cut) {
        r = sqrt(rsq);
        sigma = par[2]/r;
        f12 = 0.5*( MIN(MAX(phi[i][1],-1.0),0.0) + MIN(MAX(phi[j][1],-1.0),0.0) );
        eng_one = 4.0*par[0]*(1-par[1]*f12)*(pow(sigma,16.0) - sigma*sigma);
        eng += eng_one;
//        printf("%d %d %f\n", i+1, j+1, eng_one);
      }
    }
  }

  return eng;
}

double bend_pot(double *par, Atoms &atoms, PairList &pair)
{
  int ii,i,j;
  double r,rinv,rsq,sigma,rsq_cut,f12,uij,r12[3],ru12[3];
  double eng_one, eng = 0.0;
  double normi[3], normj[3], nti[3], ntj[3];
  double gamma_r, thetai, thetaj, alphai, alphaj;

  rsq_cut = par[4]*par[4];

  int n = pair.npair;
  int *list = pair.list;
  int *type = atoms.type;
  double **x = atoms.x;
  double **quat = atoms.quat;
  double **phi = atoms.phi;

  for (ii=0;ii<n;++ii) {
    i = list[2*ii];
    j = list[2*ii+1];
    if (type[i]==1 && type[j]==1) {
      rsq = delta(x[i], x[j], r12);
      if (rsq<rsq_cut) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        sigma = par[2]*rinv;
        f12 = 0.5*( (phi[i][1]<=0 ? -MAX(phi[i][1],-1.0) : 0.0) + (phi[j][1]<=0 ? -MAX(phi[j][1],-1.0) : 0.0) );

        ru12[0] = r12[0]*rinv;
        ru12[1] = r12[1]*rinv;
        ru12[2] = r12[2]*rinv;

        normi[0] = 2.0*(quat[i][1]*quat[i][3] + quat[i][0]*quat[i][2]);
        normi[1] = 2.0*(quat[i][2]*quat[i][3] - quat[i][0]*quat[i][1]);
        normi[2] = quat[i][0]*quat[i][0] - quat[i][1]*quat[i][1] - quat[i][2]*quat[i][2] + quat[i][3]*quat[i][3];

        normj[0] = 2.0*(quat[j][1]*quat[j][3] + quat[j][0]*quat[j][2]);
        normj[1] = 2.0*(quat[j][2]*quat[j][3] - quat[j][0]*quat[j][1]);
        normj[2] = quat[j][0]*quat[j][0] - quat[j][1]*quat[j][1] - quat[j][2]*quat[j][2] + quat[j][3]*quat[j][3];

        nti[0] = 2.0*(quat[i][1]*quat[i][2] - quat[i][0]*quat[i][3]);
        nti[1] = quat[i][0]*quat[i][0] - quat[i][1]*quat[i][1] + quat[i][2]*quat[i][2] - quat[i][3]*quat[i][3];
        nti[2] = 2.0*(quat[i][2]*quat[i][3] + quat[i][0]*quat[i][1]);

        ntj[0] = 2.0*(quat[j][1]*quat[j][2] - quat[j][0]*quat[j][3]);
        ntj[1] = quat[j][0]*quat[j][0] - quat[j][1]*quat[j][1] + quat[j][2]*quat[j][2] - quat[j][3]*quat[j][3];
        ntj[2] = 2.0*(quat[j][2]*quat[j][3] + quat[j][0]*quat[j][1]);

        thetai = dot(normi, ru12);
        thetaj = dot(normj, ru12);

        alphai = dot(nti, ru12);
        alphaj = dot(ntj, ru12);

        gamma_r = 0.25*r*par[3]*f12*(alphai*alphai + alphaj*alphaj);

        uij = (thetai - gamma_r)*(thetai - gamma_r) + (thetaj + gamma_r)*(thetaj + gamma_r);

        eng_one = 4.0*par[0]*sigma*sigma*(1+par[1]*f12)*uij;
        eng += eng_one;
//       printf("%d %d %f\n", i+1, j+1, eng_one);
      }
    }
  }

  return eng;
}

double olig_pot(double *par, Atoms &atoms, PairList &pair)
{
  int ii,i,j;
  double rsq,ainv,rsq_cut;
  double eng_one, eng = 0.0;
  double nti[3], ntj[3];

  rsq_cut = par[2]*par[2];
  ainv = 1.0/par[1];

  int n = pair.npair;
  int *list = pair.list;
  int *type = atoms.type;
  double **x = atoms.x;
  double **quat = atoms.quat;
  double **phi = atoms.phi;

  for (ii=0;ii<n;++ii) {
    i = list[2*ii];
    j = list[2*ii+1];
    if (type[i]==1 && type[j]==1) {
      rsq = dist_sq(x[i], x[j]);
      if (rsq<rsq_cut && phi[i][1]<-0.8 && phi[j][1]<-0.8) {
        nti[0] = 2.0*(quat[i][1]*quat[i][2] - quat[i][0]*quat[i][3]);
        nti[1] = quat[i][0]*quat[i][0] - quat[i][1]*quat[i][1] + quat[i][2]*quat[i][2] - quat[i][3]*quat[i][3];
        nti[2] = 2.0*(quat[i][2]*quat[i][3] + quat[i][0]*quat[i][1]);

        ntj[0] = 2.0*(quat[j][1]*quat[j][2] - quat[j][0]*quat[j][3]);
        ntj[1] = quat[j][0]*quat[j][0] - quat[j][1]*quat[j][1] + quat[j][2]*quat[j][2] - quat[j][3]*quat[j][3];
        ntj[2] = 2.0*(quat[j][2]*quat[j][3] + quat[j][0]*quat[j][1]);

        eng_one = -par[0]*dot(nti,ntj);
        if (par[1]!=1.0) eng_one *= ainv*pow(dot(nti,ntj),par[1]-1.0);
        eng += eng_one;
      }
    }
  }

  return eng;
}

double cc_pot(double *par, Atoms &atoms, PairList &pair=PVAL)
{
  int i;
  double zeta_m, zeta_b;
  double eng_one, eng = 0.0;

  int n = atoms.n_mem;
  double **phi = atoms.phi;

  for (i=0;i<n;++i) {
    zeta_m = 0.5*(phi[i][0]-1.0) - par[1];
    zeta_b = 0.5*(phi[i][1]-1.0) - par[1];
    eng_one = -par[0]*zeta_m*zeta_b;
    eng += eng_one;
  }

  return eng;
}

double ic_pot(double *par, Atoms &atoms, PairList &pair)
{
  int ii,i,j;
  double r,rinv,rsq,rsq_cut,f12,uij,r12[3],ru12[3];
  double eng_one, eng = 0.0;
  double normi[3], normj[3], nti[3], ntj[3];
  double gamma_r, thetai, thetaj, alphai, alphaj;

  rsq_cut = par[3]*par[3];

  int nlist = pair.npair;
  int *list = pair.list;
  int n_mem = atoms.n_mem;
  int *type = atoms.type;
  double **x = atoms.x;
  double **quat = atoms.quat;
  double **phi = atoms.phi;

  for (i=0;i<n_mem;++i) {
    eng_one = par[0]*MIN(MAX(phi[i][1],-1.0),1.0);
    eng += eng_one;
//    printf("%d %f %f\n", i+1, phi[i][1], eng_one);
  } 

  for (ii=0;ii<nlist;++ii) {
    i = list[2*ii];
    j = list[2*ii+1];
    if (type[i]==1 && type[j]==1) {
      rsq = delta(x[i], x[j], r12);
      if (rsq<rsq_cut) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        f12 = 0.5*( (phi[i][1]<=0 ? -MAX(phi[i][1],-1.0) : 0.0) + (phi[j][1]<=0 ? -MAX(phi[j][1],-1.0) : 0.0) );

        ru12[0] = r12[0]*rinv;
        ru12[1] = r12[1]*rinv;
        ru12[2] = r12[2]*rinv;

        normi[0] = 2.0*(quat[i][1]*quat[i][3] + quat[i][0]*quat[i][2]);
        normi[1] = 2.0*(quat[i][2]*quat[i][3] - quat[i][0]*quat[i][1]);
        normi[2] = quat[i][0]*quat[i][0] - quat[i][1]*quat[i][1] - quat[i][2]*quat[i][2] + quat[i][3]*quat[i][3];

        normj[0] = 2.0*(quat[j][1]*quat[j][3] + quat[j][0]*quat[j][2]);
        normj[1] = 2.0*(quat[j][2]*quat[j][3] - quat[j][0]*quat[j][1]);
        normj[2] = quat[j][0]*quat[j][0] - quat[j][1]*quat[j][1] - quat[j][2]*quat[j][2] + quat[j][3]*quat[j][3];

        nti[0] = 2.0*(quat[i][1]*quat[i][2] - quat[i][0]*quat[i][3]);
        nti[1] = quat[i][0]*quat[i][0] - quat[i][1]*quat[i][1] + quat[i][2]*quat[i][2] - quat[i][3]*quat[i][3];
        nti[2] = 2.0*(quat[i][2]*quat[i][3] + quat[i][0]*quat[i][1]);

        ntj[0] = 2.0*(quat[j][1]*quat[j][2] - quat[j][0]*quat[j][3]);
        ntj[1] = quat[j][0]*quat[j][0] - quat[j][1]*quat[j][1] + quat[j][2]*quat[j][2] - quat[j][3]*quat[j][3];
        ntj[2] = 2.0*(quat[j][2]*quat[j][3] + quat[j][0]*quat[j][1]);

        thetai = dot(normi, ru12);
        thetaj = dot(normj, ru12);

        alphai = dot(nti, ru12);
        alphaj = dot(ntj, ru12);

        gamma_r = 0.25*r*par[2]*f12*(alphai*alphai + alphaj*alphaj);

        uij = (thetai - gamma_r)*(thetai - gamma_r) + (thetaj + gamma_r)*(thetaj + gamma_r);

        eng_one = -par[1]*(MIN(MAX(phi[i][1],-1.0),1.0) + MIN(MAX(phi[j][1],-1.0),1.0))*uij;
        eng += eng_one;
//        printf("%d %d %f\n", i+1, j+1, eng_one);
      }
    }
  }

  return eng;
}

void calc_rho(double rc, Atoms &atoms, PairList &pair, double **rho)
{
  int ii,i,j;
  double rsq,r,rcsq,rcinv,sigma,w;

  rcsq = rc*rc;
  rcinv = 1.0/rc;

  int nlist = pair.npair;
  int *list = pair.list;
  int n_atom = atoms.n_atom;
  int *type = atoms.type;
  double **x = atoms.x;

  for (i=0;i<n_atom;i++) {
    if (type[i]==1) rho[i][0] = 1.0;
    else rho[i][0] = 0.0;
    rho[i][1] = 1.0;
  }

  for (ii=0;ii<nlist;++ii) {
    i = list[2*ii];
    j = list[2*ii+1];
    rsq = dist_sq(x[i], x[j]);
    if (rsq<rcsq) {
      r = sqrt(rsq);
      sigma = r*rcinv;
      w = (1-sigma)*(1-sigma)*(1-sigma)*(1+3*sigma);
      if (type[i]==1 && type[j]==1) {
        rho[i][0] += w;
        rho[j][0] += w;
      }
      rho[i][1] += w;
      rho[j][1] += w;
    }
  }
}

void calc_phi_gard(double rc, Atoms &atoms, double **rho, PairList &pair, double **xi_m, double **xi_b)
{
  int ii,i,j;
  double r12[3];
  double rsq,r,rcsq,rcinv,sigma,dw,rho_ave,xi;

  rcsq = rc*rc;
  rcinv = 1.0/rc;

  int nlist = pair.npair;
  int *list = pair.list;
  int n_atom = atoms.n_atom;
  int *type = atoms.type;
  double **x = atoms.x;
  double **phi = atoms.phi;

  for (i=0;i<n_atom;i++) {
    xi_m[i][0] = 0.0;
    xi_m[i][1] = 0.0;
    xi_m[i][2] = 0.0;

    xi_b[i][0] = 0.0;
    xi_b[i][1] = 0.0;
    xi_b[i][2] = 0.0;
  }

  for (ii=0;ii<nlist;++ii) {
    i = list[2*ii];
    j = list[2*ii+1];
    rsq = delta(x[i], x[j], r12);
    if (rsq<rcsq) {
      r = sqrt(rsq);
      sigma = r*rcinv;
      dw = -12*(1-sigma)*(1-sigma)*rcinv*rcinv;
      if (type[i]==1 && type[j]==1) {
        rho_ave = 0.5*(rho[i][0] + rho[j][0]);
        xi = (phi[j][0] - phi[i][0])*dw/rho_ave;

        xi_m[i][0] += xi*r12[0];
        xi_m[i][1] += xi*r12[1];
        xi_m[i][2] += xi*r12[2];

        xi_m[j][0] += xi*r12[0];
        xi_m[j][1] += xi*r12[1];
        xi_m[j][2] += xi*r12[2];
      }
      rho_ave = 0.5*(rho[i][1] + rho[j][1]);
      xi = (phi[j][1] - phi[i][1])*dw/rho_ave;

      xi_b[i][0] += xi*r12[0];
      xi_b[i][1] += xi*r12[1];
      xi_b[i][2] += xi*r12[2];

      xi_b[j][0] += xi*r12[0];
      xi_b[j][1] += xi*r12[1];
      xi_b[j][2] += xi*r12[2];
    }
  }
}

void calc_phi_gard_sq(Atoms &atoms, double **xi_m, double **xi_b, double **xi_sq)
{
  int i;

  int *type = atoms.type;
  int n_atom = atoms.n_atom;

  for (i=0;i<n_atom;i++) {
    if (type[i]==1) xi_sq[i][0] = dot(xi_m[i], xi_m[i]);
    else xi_sq[i][0] = 0.0;
    xi_sq[i][1] = dot(xi_b[i], xi_b[i]);
  }
}

double grad_pots(double *par, Atoms &atoms, PairList &pair)
{
  int i;
  double eng_m, eng_b;
  double rcut;

  rcut = par[2];

  int n_atom = atoms.n_atom;
  int *type = atoms.type;

  double **rho = (double **)malloc(n_atom * sizeof(double *));
  double **xi_m = (double **)malloc(n_atom * sizeof(double *));
  double **xi_b = (double **)malloc(n_atom * sizeof(double *));
  for (i = 0; i < n_atom; i++) {
     rho[i] = (double *)malloc(2 * sizeof(double));
     xi_m[i] = (double *)malloc(3 * sizeof(double));
     xi_b[i] = (double *)malloc(3 * sizeof(double));
  }

  calc_rho(rcut, atoms, pair, rho);
  calc_phi_gard(rcut, atoms, rho, pair, xi_m, xi_b);

  eng_m = eng_b = 0.0;
  for (i = 0; i < n_atom; i++) {
    if (type[i]==1) eng_m += dot(xi_m[i], xi_m[i]);
    eng_b += dot(xi_b[i], xi_b[i]);
  }
  
  eng_m *= par[0];
  eng_b *= par[1];

  printf("eng_m = %f eng_b = %f\n", eng_m, eng_b);

  for (i = 0; i < n_atom; i++) {
    free(rho[i]);
    free(xi_m[i]);
    free(xi_b[i]);
  }
  free(rho);
  free(xi_m);
  free(xi_b);

  return eng_m + eng_b;
}

double grad_pots(double *par, Atoms &atoms, PairList &pair, double **rho, double **xi_m, double **xi_b)
{
  int i;
  double eng_m, eng_b;
  double rcut;

  rcut = par[2];

  int n_atom = atoms.n_atom;
  int *type = atoms.type;

  eng_m = eng_b = 0.0;
  for (i = 0; i < n_atom; i++) {
    if (type[i]==1) eng_m += dot(xi_m[i], xi_m[i]);
    eng_b += dot(xi_b[i], xi_b[i]);
  }

  eng_m *= par[0];
  eng_b *= par[1];

  return eng_m + eng_b;
}

double mem_comp_pot(double *par, Atoms &atoms, PairList &pair=PVAL)
{
  int i;
  double eng_one, eng = 0.0;

  int n = atoms.n_mem;
  int *type = atoms.type;
  double **phi = atoms.phi;

  for (i=0;i<n;i++) {
    if (type[i]==1) {
      eng_one = 0.1*par[0]*pow(phi[i][0],10.0);
      eng += eng_one;
    }
  }

  return eng;
}

double prot_comp_pot(double *par, Atoms &atoms, PairList &pair=PVAL)
{
  int i;
  double eng_one, eng = 0.0;
  double phi_sq;

  int n = atoms.n_atom;
  double **phi = atoms.phi;

  for (i=0;i<n;i++) {
    phi_sq = phi[i][1]*phi[i][1];
    eng_one = par[0]*phi_sq*(phi_sq*phi_sq + 1.0);
    eng += eng_one;
  }

  return eng;
}

typedef double (*function_type1)(double*, Atoms&);
typedef double (*function_type2)(double*, Atoms&, PairList&);

void calc_force_one(int i, double dx, double *par, function_type2 pot, Atoms &atoms, PairList &pair, double *f=NULL, int pflag=1)
{
  int j;
  double xtmp[3], eng1, eng2, fone[3];

  double **x = atoms.x;
  
  xtmp[0] = x[i][0];
  xtmp[1] = x[i][1];
  xtmp[2] = x[i][2];

  for (j=0;j<3;j++) {
    x[i][j] = xtmp[j] - dx;
    eng1 = pot(par, atoms, pair);
    x[i][j] = xtmp[j] + dx;
    eng2 = pot(par, atoms, pair);
    x[i][j] = xtmp[j];
    fone[j] = -0.5*(eng2 - eng1)/dx;
  }

  if (pflag) printf("f[%d]={%.8f, %.8f, %.8f}\n", i+1, fone[0], fone[1], fone[2]);

  if (f!=NULL) {
    f[0] = fone[0];
    f[1] = fone[1];
    f[2] = fone[2];
  }
}

void calc_force(double dx, double *par, function_type2 pot, Atoms &atoms, PairList &pair, double **f=NULL, int pflag=1)
{
  int i,j;
  double xtmp[3], eng1, eng2, fone[3];

  int n = pair.natom;
  double **x = atoms.x;

  for (i=0;i<n;i++) {
  
    xtmp[0] = x[i][0];
    xtmp[1] = x[i][1];
    xtmp[2] = x[i][2];

    for (j=0;j<3;j++) {
      x[i][j] = xtmp[j] - dx;
      eng1 = pot(par, atoms, pair);
      x[i][j] = xtmp[j] + dx;
      eng2 = pot(par, atoms, pair);
      x[i][j] = xtmp[j];
      fone[j] = -0.5*(eng2 - eng1)/dx;
    }

    if (pflag) printf("%d %.8f %.8f %.8f\n", i+1, fone[0], fone[1], fone[2]);

    if (f!=NULL) {
      f[i][0] = fone[0];
      f[i][1] = fone[1];
      f[i][2] = fone[2];
    }
  }
}

inline void dquat(double *q, double *q0, double sina, double cosa, int ix)
{
  q[0] = cosa*q0[0];
  q[1] = cosa*q0[1];
  q[2] = cosa*q0[2];
  q[3] = cosa*q0[3];
  switch (ix) {
    case 0:
      q[0] -= sina*q0[1];
      q[1] += sina*q0[0];
      q[2] -= sina*q0[3];
      q[3] += sina*q0[2];
      break;
    case 1:
      q[0] -= sina*q0[2];
      q[1] += sina*q0[3];
      q[2] += sina*q0[0];
      q[3] -= sina*q0[1];
      break;
    case 2:
      q[0] -= sina*q0[3];
      q[1] -= sina*q0[2];
      q[2] += sina*q0[1];
      q[3] += sina*q0[0];
      break;
    default:
      printf("Error in dquat function\n");
      exit(0);
  }
}

void calc_torque_one(int i, double dp, double *par, function_type2 pot, Atoms &atoms, PairList &pair, double *t=NULL, int pflag=1)
{
  int j;
  double qtmp[4], eng1, eng2, tone[3];
  double sina, cosa;

  double **quat  = atoms.quat;

  sina = sin(dp/2.0);
  cosa = cos(dp/2.0);
  
  qtmp[0] = quat[i][0];
  qtmp[1] = quat[i][1];
  qtmp[2] = quat[i][2];
  qtmp[3] = quat[i][3];

  for (j=0;j<3;j++) {
    dquat(quat[i], qtmp, -sina, cosa, j);
    eng1 = pot(par, atoms, pair);
    dquat(quat[i], qtmp, sina, cosa, j);
    eng2 = pot(par, atoms, pair);
    tone[j] = -0.5*(eng2 - eng1)/dp;
  }

  quat[i][0] = qtmp[0];
  quat[i][1] = qtmp[1];
  quat[i][2] = qtmp[2];
  quat[i][3] = qtmp[3];

  if (pflag) printf("t[%d]={%.8f, %.8f, %.8f}\n", i+1, tone[0], tone[1], tone[2]);

  if (t!=NULL) {
    t[0] = tone[0];
    t[1] = tone[1];
    t[2] = tone[2];
  }
}

void calc_torque(double dp, double *par, function_type2 pot, Atoms &atoms, PairList &pair, double **t=NULL, int pflag=1)
{
  int i,j;
  double qtmp[4], eng1, eng2, tone[3];
  double sina, cosa;

  int n = pair.natom;
  double **quat  = atoms.quat;

  sina = sin(dp/2.0);
  cosa = cos(dp/2.0);

  for (i=0;i<n;i++) {
    qtmp[0] = quat[i][0];
    qtmp[1] = quat[i][1];
    qtmp[2] = quat[i][2];
    qtmp[3] = quat[i][3];

    for (j=0;j<3;j++) {
      dquat(quat[i], qtmp, -sina, cosa, j);
      eng1 = pot(par, atoms, pair);
      dquat(quat[i], qtmp, sina, cosa, j);
      eng2 = pot(par, atoms, pair);
      tone[j] = -0.5*(eng2 - eng1)/dp;
    }

    quat[i][0] = qtmp[0];
    quat[i][1] = qtmp[1];
    quat[i][2] = qtmp[2];
    quat[i][3] = qtmp[3];

    if (pflag) printf("%d %.8f %.8f %.8f\n", i+1, tone[0], tone[1], tone[2]);

    if (t!=NULL) {
      t[i][0] = tone[0];
      t[i][1] = tone[1];
      t[i][2] = tone[2];
    }
  }
}

void calc_dphi_one(int i, double dp, double *par, function_type2 pot, Atoms &atoms, PairList &pair, double *dphi=NULL, int pflag=1)
{
  int j;
  double phi_tmp[3], eng1, eng2, dphi_one[3];
  double factor = 0.5*phi_g/dp;

  double **phi = atoms.phi;
  
  phi_tmp[0] = phi[i][0];
  phi_tmp[1] = phi[i][1];

  for (j=0;j<2;j++) {
    phi[i][j] = phi_tmp[j] - dp;
    eng1 = pot(par, atoms, pair);
    phi[i][j] = phi_tmp[j] + dp;
    eng2 = pot(par, atoms, pair);
    phi[i][j] = phi_tmp[j];
    dphi_one[j] = -factor*(eng2 - eng1);
  }

  if (pflag) printf("dphi[%d]={%.8f, %.8f}\n", i+1, dphi_one[0], dphi_one[1]);

  if (dphi!=NULL) {
    dphi[0] = dphi_one[0];
    dphi[1] = dphi_one[1];
  }
}

void calc_dphi(double dp, double *par, function_type2 pot, Atoms &atoms, PairList &pair, double **dphi=NULL, int pflag=1)
{
  int i,j;
  double phi_tmp[3], eng1, eng2, dphi_one[3];
  double factor = 0.5*phi_g/dp;
  
  int n = pair.natom;
  double **phi = atoms.phi;
  
  for (i=0;i<n;i++) {
    phi_tmp[0] = phi[i][0];
    phi_tmp[1] = phi[i][1];

    for (j=0;j<2;j++) {
      phi[i][j] = phi_tmp[j] - dp;
      eng1 = pot(par, atoms, pair);
      phi[i][j] = phi_tmp[j] + dp;
      eng2 = pot(par, atoms, pair);
      phi[i][j] = phi_tmp[j];
      dphi_one[j] = -factor*(eng2 - eng1);
    }

    if (pflag) printf("%d %.8f %.8f\n", i+1, dphi_one[0], dphi_one[1]);

    if (dphi!=NULL) {
      dphi[i][0] = dphi_one[0];
      dphi[i][1] = dphi_one[1];
    }
  }
}

void calc_dphi_grad_pots(double dp, double *par, function_type2 pot, Atoms &atoms, PairList &pair, double **dphi=NULL, int pflag=1)
{
  int i,j;
  double phi_tmp[3], eng1, eng2, dphi_one[3];
  double factor = 0.5*phi_g/dp;
  double rcut = par[2];

  int n = pair.natom;
  double **phi = atoms.phi;

  double **rho = (double **)malloc(n * sizeof(double *));
  double **xi_m = (double **)malloc(n * sizeof(double *));
  double **xi_b = (double **)malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
     rho[i] = (double *)malloc(2 * sizeof(double));
     xi_m[i] = (double *)malloc(3 * sizeof(double));
     xi_b[i] = (double *)malloc(3 * sizeof(double));
  }

  calc_rho(rcut, atoms, pair, rho);

  for (i=0;i<n;i++) {
    phi_tmp[0] = phi[i][0];
    phi_tmp[1] = phi[i][1];

    for (j=0;j<2;j++) {
      phi[i][j] = phi_tmp[j] - dp;
      calc_phi_gard(rcut, atoms, rho, pair, xi_m, xi_b);
      eng1 = grad_pots(par, atoms, pair, rho, xi_m, xi_b);
      phi[i][j] = phi_tmp[j] + dp;
      calc_phi_gard(rcut, atoms, rho, pair, xi_m, xi_b);
      eng2 = grad_pots(par, atoms, pair, rho, xi_m, xi_b);
      phi[i][j] = phi_tmp[j];
      dphi_one[j] = -factor*(eng2 - eng1);
    }

    if (pflag) printf("%d %.8f %.8f\n", i+1, dphi_one[0], dphi_one[1]);

    if (dphi!=NULL) {
      dphi[i][0] = dphi_one[0];
      dphi[i][1] = dphi_one[1];
    }
  }

  for (i = 0; i < n; i++) {
    free(rho[i]);
    free(xi_m[i]);
    free(xi_b[i]);
  }
  free(rho);
  free(xi_m);
  free(xi_b);
}

int main()
{
  int i,n;
  int n_mem = 5882;
  int n_prot = 89134;
  int n_tot = n_mem + n_prot;

  double lucy_pars[] = {1.195, 0.239, 136.0};
  double lj216_pars[] = {1.553537, 0.5, 68.0, 136.0};
  double bend_pars[] = {2.50956, 0.5, 68.0, 0.009, 136.0};
  double olig_pars[] = {0.5, 1.0, 136.0};
  double ic_pars[] = {0.06, 1.195, 0.009, 136.0};
  double cc_pars[] = {9.56022944, 2.0};
  double grad_pars[] = {0.239, 0.239, 136.0};
  double mem_comp_pars[] = {0.000239};
  double prot_comp_pars[] = {0.000239};

  char xyz_file[] = "xyz";

  FILE *file = fopen(xyz_file, "r");

  Atoms atoms(n_mem, n_tot);

  n = 0;
  int id, ty, nf;
  double q[4], p[2], xyz[3];
  while (fscanf(file,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",&id,&ty,&q[0],&q[1],&q[2],&q[3],&p[0],&p[1],&xyz[0],&xyz[1],&xyz[2]) != EOF) {

    if (n>=n_tot) {
      printf("n: %d nf: %d\n",n,nf);
      printf("XYZ file format error!\n");
      return -1;
    }

    atoms.type[n] = ty;
    atoms.quat[n][0] = q[0];
    atoms.quat[n][1] = q[1];
    atoms.quat[n][2] = q[2];
    atoms.quat[n][3] = q[3];
    atoms.phi[n][0] = p[0];
    atoms.phi[n][1] = p[1];
    atoms.x[n][0] = xyz[0];
    atoms.x[n][1] = xyz[1];
    atoms.x[n][2] = xyz[2];

    n++;
  }

  if (n!=n_tot) {
    printf("END n: %d n_tot: %d\n",n,n_tot);
    printf("XYZ file format error!\n");
    return -1;
  }

  fclose(file);

  double dx = 0.00001;
  double dp = 0.00001;
  double dph = 0.00001;

  // Making pair lists

  PairList plist_mem;
  PairList plist_all;
  double rcut = 136.0;
  double rcutsq = (rcut+2.0*dx)*(rcut+2.0*dx);
  printf("Makeing mem pair list ...\n");
  int npair_mem = plist_mem.make_pair_list(rcutsq, atoms, n_mem);
  printf("Makeing all pair list ...\n");
//  int npair_all = plist_all.make_pair_list(rcutsq, atoms, n_tot);

//  printf("npair_mem=%d npair_all=%d\n",npair_mem,npair_all);

  double eng=0.0;
  eng = lj216_pot(lj216_pars, atoms, plist_mem);
  printf("LJ216 Energy: %f\n", eng);

  eng = lucy_pot(lucy_pars, atoms, plist_all);
  printf("Lucy Energy: %f\n", eng);

  eng = bend_pot(bend_pars, atoms, plist_mem);
  printf("Bending Energy: %f\n", eng);

  eng = olig_pot(olig_pars, atoms, plist_mem);
  printf("Oligomerization Energy: %f\n", eng);

  eng = ic_pot(ic_pars, atoms, plist_mem);
  printf("Intrinsic Curvature Energy: %f\n", eng);

  eng = cc_pot(cc_pars, atoms);
  printf("Composition Coupling Energy: %f\n", eng);

  eng = grad_pots(grad_pars, atoms, plist_all);
  printf("Gradient Energy: %f\n", eng);

  eng = mem_comp_pot(mem_comp_pars, atoms);
  printf("Memory Composition Energy: %f\n", eng);

  eng = prot_comp_pot(prot_comp_pars, atoms);
  printf("Protenin Composition Energy: %f\n", eng);

//  calc_force_one(0, dx, bend_pars, &bend_pot, atoms, plist_mem);
//  calc_force(dx, bend_pars, &bend_pot, atoms, plist_mem);
//  calc_torque_one(0, dp, bend_pars, &bend_pot, atoms, plist_mem);
//  calc_torque(dp, bend_pars, &bend_pot, atoms, plist_mem);
//  calc_dphi_one(0, dph, bend_pars, &bend_pot, atoms, plist_mem);
//  calc_dphi(dph, bend_pars, &bend_pot, atoms, plist_mem);

//  calc_force(dx, lj216_pars, &lj216_pot, atoms, plist_mem);
  calc_dphi(dph, lj216_pars, &lj216_pot, atoms, plist_mem);

//  calc_force_one(0, dx, ic_pars, &ic_pot, atoms, plist_mem);
//  calc_torque_one(0, dp, ic_pars, &ic_pot, atoms, plist_mem);
//  calc_dphi_one(0, dph, ic_pars, &ic_pot, atoms, plist_mem);

//  calc_torque(dp, olig_pars, &olig_pot, atoms, plist_mem);

//  calc_dphi_one(0, dph, cc_pars, &cc_pot, atoms, plist_mem);
//  calc_dphi(dph, cc_pars, &cc_pot, atoms, plist_mem);

//  calc_dphi(dph, mem_comp_pars, &mem_comp_pot, atoms, plist_mem);
//  calc_dphi(dph, prot_comp_pars, &prot_comp_pot, atoms, plist_all);

//  calc_dphi_grad_pots(dph, grad_pars, &grad_pots, atoms, plist_all);
//  calc_dphi_grad_pots(dph, grad_pars, &grad_pots, atoms, plist_all);


/*  double **rho = (double **)malloc(n_tot * sizeof(double *));
  double **xi_m = (double **)malloc(n_tot * sizeof(double *));
  double **xi_b = (double **)malloc(n_tot * sizeof(double *));
  for (i = 0; i < n_tot; i++) {
     rho[i] = (double *)malloc(2 * sizeof(double));
     xi_m[i] = (double *)malloc(3 * sizeof(double));
     xi_b[i] = (double *)malloc(3 * sizeof(double));
  }

//  calc_rho(rcut, atoms, plist_all, rho);
//  calc_phi_gard(rcut, atoms, rho, plist_all, xi_m, xi_b);
//  calc_rho(rcut, type, x, plist_all, n_tot, npair_all, rho);
//  calc_phi_gard(rcut, type, x, phi, rho, plist_all, n_tot, npair_all, xi_m, xi_b);

//  for (i = 0; i < n_tot; i++) {
//    printf("%d %f %f\n", i+1, rho[i][0], rho[i][1]);
//    printf("%d %f %f %f %f %f %f\n", i+1, xi_m[i][0], xi_m[i][1], xi_m[i][2], xi_b[i][0], xi_b[i][1], xi_b[i][2]);
//  }


  for (i = 0; i < n_tot; i++) {
    free(rho[i]);
    free(xi_m[i]);
    free(xi_b[i]);
  }
  free(rho);
  free(xi_m);
  free(xi_b);*/

  return 0;
}
