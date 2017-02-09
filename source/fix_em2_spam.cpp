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
#include "fix_em2_spam.h"
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

FixEM2SPAM::FixEM2SPAM(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg) 
{
  if (narg != 4) error->all(FLERR,"Illegal fix em2/spam command");

  allocated = 0;

  int len = strlen(arg[3]) + 1;
  conf_file = new char[len];
  strcpy(conf_file, arg[3]);
}

/* ---------------------------------------------------------------------- */

FixEM2SPAM::~FixEM2SPAM()
{
  if (allocated) {
    memory->destroy(inertia);
    memory->destroy(quat_flag);
    memory->destroy(mem_flag);
    memory->destroy(prot_flag);
  }
}

/* ---------------------------------------------------------------------- */

void FixEM2SPAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(inertia,n+1,3,"fix:inertia");
  memory->create(quat_flag,n+1,"fix:quat_flag");
  memory->create(mem_flag,n+1,"fix:mem_flag");
  memory->create(prot_flag,n+1,"fix:protmem_flag");
}

void FixEM2SPAM::init()
{
  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec)
    error->all(FLERR,"Compute em2/spam requires atom style em2");

  FixNVE::init();

  dtv2 = 0.5 * dtv;

  if (!allocated) allocate();
  read_conf_file(conf_file);
}

/* ---------------------------------------------------------------------- */

void FixEM2SPAM::initial_integrate(int vflag)
{
  int itype;

  double **phi = avec->phi;
  double **phi_half = avec->phi_half;
  double **dphi = avec->dphi;
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

void FixEM2SPAM::final_integrate()
{
  int itype;

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **phi_half = avec->phi_half;
  double **dphi = avec->dphi;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      itype = type[i];

      // Final SPAM integration

      if (mem_flag[itype]) phi_half[i][0] += dtv2*dphi[i][0];
      if (prot_flag[itype]) phi_half[i][1] += dtv2*dphi[i][1];
    }
}

/* ---------------------------------------------------------------------- */

void FixEM2SPAM::reset_dt()
{
  FixNVE::reset_dt();

  dtv2 = 0.5 * dtv;
}

/* ---------------------------------------------------------------------- */

void FixEM2SPAM::read_conf_file(char *filename)
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
  if (!file) error->all(FLERR,"Fix em2/spam: Error opening configuration file!");

  if (me==0) print_log("Reading integrator configuration file ...\n");

  while ( fgets ( ln, sizeof ln, file ) != NULL ) {
    line = trim(ln);

    if (line[0]=='#' || isEmptyString(line)) continue;

    // Split the parameter line into strings
    narg = 0;
    arg[narg] = strtok(line, " \t\n");
    while ( arg[narg]!=NULL ) {
      narg++;
      if (narg>=16) error->all(FLERR,"Fix em2/spam: Configuration file line is too long!");
      arg[narg] = strtok(NULL, " \t\n");
    }

    if (narg!=7) error->all(FLERR,"Fix em2/spam: Configuration file format error");

    I[0] = atof(arg[1]);
    I[1] = atof(arg[2]);
    I[2] = atof(arg[3]);
    flag1 = atoi(arg[4]);
    flag2 = atoi(arg[5]);
    flag3 = atoi(arg[6]);

    force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);

    if ((flag1!=0 && flag1!=1) || (flag2!=0 && flag2!=1) || (flag3!=0 && flag3!=1))
      error->all(FLERR,"Fix em2/spam: Configuration file format error");

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

    if (count==0) error->all(FLERR,"Fix em2/spam: Configuration file format error");
  }

  fclose(file);

  for (i=1;i<=n;i++) {
    if (!type_flag) error->all(FLERR,"Fix em2/spam: Parameters for all types must be set!");
  }

  memory->destroy(type_flag);
}

inline void FixEM2SPAM::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

char *FixEM2SPAM::ltrim(char *s)
{     
  while(isspace(*s)) s++;     
  return s; 
}  

char *FixEM2SPAM::rtrim(char *s)
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

char *FixEM2SPAM::trim(char *s)
{     
  return rtrim(ltrim(s));  
}

bool FixEM2SPAM::isEmptyString(char *str)
{
  int len = strlen(str);
  
  if (len==0) return true;
  
  for (int i=0;i<len;++i) {
    if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') return false;
  }

  return true;
}
