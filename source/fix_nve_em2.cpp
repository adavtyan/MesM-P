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
#include "fix_nve_em2.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_em2.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEEM2::FixNVEEM2(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg) 
{
  if (narg != 4) error->all(FLERR,"Illegal fix nve/em2 command");

  allocated = 0;

  int len = strlen(arg[3]) + 1;
  conf_file = new char[len];
  strcpy(conf_file, arg[3]);
}

/* ---------------------------------------------------------------------- */

FixNVEEM2::~FixNVEEM2()
{
  if (allocated) {
    memory->destroy(inertia);
    memory->destroy(quat_flag);
    memory->destroy(mem_flag);
    memory->destroy(prot_flag);
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEEM2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(inertia,n+1,3,"fix:inertia");
  memory->create(quat_flag,n+1,"fix:quat_flag");
  memory->create(mem_flag,n+1,"fix:mem_flag");
  memory->create(prot_flag,n+1,"fix:protmem_flag");
}

void FixNVEEM2::init()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec) avec = (AtomVecEM2 *) atom->style_match("em2_angle");
  if (!avec) error->all(FLERR,"Compute nve/em2 requires atom style em2");

  FixNVE::init();

  dtv2 = 0.5 * dtv;

  if (!allocated) allocate();
  read_conf_file(conf_file);
}

/* ---------------------------------------------------------------------- */

void FixNVEEM2::initial_integrate(int vflag)
{
  double dtfm;
  double omega[3];
  double *quat,*inertia_i;
  int itype;

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **phi = avec->phi;
  double **phi_half = avec->phi_half;
  double **dphi = avec->dphi;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  // Should be solved by reset_dt()
  // dtv2 = 0.5 * dtv; 

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      itype = type[i];

      if (rmass) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[itype];
      }
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // update angular momentum by 1/2 step

      if (quat_flag[itype]) {
        angmom[i][0] += dtf * torque[i][0];
        angmom[i][1] += dtf * torque[i][1];
        angmom[i][2] += dtf * torque[i][2];

        // compute omega at 1/2 step from angmom at 1/2 step and current q
        // update quaternion a full step via Richardson iteration
        // returns new normalized quaternion

        inertia_i = inertia[itype];
        quat = bonus[i].quat;
        MathExtra::mq_to_omega(angmom[i],quat,inertia_i,omega);
        MathExtra::richardson(quat,angmom[i],omega,inertia_i,dtv2);
      }

      // Initial SPAM integration

      if (mem_flag[itype]) {
        phi[i][0] = phi_half[i][0] + dtv*dphi[i][0];
        phi_half[i][0] += dtv2*dphi[i][0];
      }

      if (prot_flag[itype]) {
        phi[i][1] = phi_half[i][1] + dtv*dphi[i][1];
        phi_half[i][1] += dtv2*dphi[i][1];
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixNVEEM2::final_integrate()
{
  double dtfm;
  int itype;

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **phi_half = avec->phi_half;
  double **dphi = avec->dphi;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      itype = type[i];

      if (rmass) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[itype];
      }
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      if (quat_flag[itype]) {
        angmom[i][0] += dtf * torque[i][0];
        angmom[i][1] += dtf * torque[i][1];
        angmom[i][2] += dtf * torque[i][2];
      }

      // Final SPAM integration

      if (mem_flag[itype]) phi_half[i][0] += dtv2*dphi[i][0];
      if (prot_flag[itype]) phi_half[i][1] += dtv2*dphi[i][1];
    }
}

/* ---------------------------------------------------------------------- */

void FixNVEEM2::reset_dt()
{
  FixNVE::reset_dt();

  dtv2 = 0.5 * dtv;
}

/* ---------------------------------------------------------------------- */

void FixNVEEM2::read_conf_file(char *filename)
{
  int i;
  int ilo,ihi,count;
  int flag1, flag2, flag3;
  double I[3];
  int *type_flag;

  int me = comm->me;

  FILE *file;
  int narg=0;
  char ln[1024], *line, *arg[16];

  int n=atom->ntypes;
  memory->create(type_flag,n+1,"fix:type_flag");

  for (i=1;i<=n;i++) type_flag[i] = 0;

  file = fopen(filename,"r");
  if (!file) error->all(FLERR,"Fix nve/em2: Error opening configuration file!");

  if (me==0) print_log("Reading integrator configuration file ...\n");

  while ( fgets ( ln, sizeof ln, file ) != NULL ) {
    line = trim(ln);

    if (line[0]=='#' || isEmptyString(line)) continue;

    // Split the parameter line into strings
    narg = 0;
    arg[narg] = strtok(line, " \t\n");
    while ( arg[narg]!=NULL ) {
      narg++;
      if (narg>=16) error->all(FLERR,"Fix nve/em2: Configuration file line is too long!");
      arg[narg] = strtok(NULL, " \t\n");
    }

    if (narg!=7) error->all(FLERR,"Fix nve/em2: Configuration file format error");

    I[0] = atof(arg[1]);
    I[1] = atof(arg[2]);
    I[2] = atof(arg[3]);
    flag1 = atoi(arg[4]);
    flag2 = atoi(arg[5]);
    flag3 = atoi(arg[6]);

    force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);

    if ((flag1!=0 && flag1!=1) || (flag2!=0 && flag2!=1) || (flag3!=0 && flag3!=1))
      error->all(FLERR,"Fix nve/em2: Configuration file format error");

    count = 0;
    for (i = ilo; i <= ihi; i++) {
      type_flag[i] = 1;
      inertia[i][0] = I[0];
      inertia[i][1] = I[1];
      inertia[i][2] = I[2];
      quat_flag[i] = flag1;
      mem_flag[i] = flag2;
      prot_flag[i] = flag3;
      count++;
    }

    if (count==0) error->all(FLERR,"Fix nve/em2: Configuration file format error");
  }

  fclose(file);

  for (i=1;i<=n;i++) {
    if (!type_flag) error->all(FLERR,"Fix nve/em2: Parameters for all types must be set!");
  }

  memory->destroy(type_flag);
}

inline void FixNVEEM2::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

char *FixNVEEM2::ltrim(char *s)
{     
  while(isspace(*s)) s++;     
  return s; 
}  

char *FixNVEEM2::rtrim(char *s)
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

char *FixNVEEM2::trim(char *s)
{     
  return rtrim(ltrim(s));  
}

bool FixNVEEM2::isEmptyString(char *str)
{
  int len = strlen(str);
  
  if (len==0) return true;
  
  for (int i=0;i<len;++i) {
    if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') return false;
  }

  return true;
}
