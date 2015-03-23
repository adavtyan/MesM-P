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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "compute_temp_em2.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_em2.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{ROTATE,ALL};

/* ---------------------------------------------------------------------- */

ComputeTempEM2::ComputeTempEM2(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute temp/em2 command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  allocated = 0;

  int len = strlen(arg[3]) + 1;
  conf_file = new char[len];
  strcpy(conf_file, arg[3]);

  tempbias = 0;
  id_bias = NULL;
  mode = ALL;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg) 
        error->all(FLERR,"Illegal compute temp/em2 command");
      tempbias = 1;
      int n = strlen(arg[iarg+1]) + 1;
      id_bias = new char[n];
      strcpy(id_bias,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dof") == 0) {
      if (iarg+2 > narg) 
        error->all(FLERR,"Illegal compute temp/em2 command");
      if (strcmp(arg[iarg+1],"rotate") == 0) mode = ROTATE;
      else if (strcmp(arg[iarg+1],"all") == 0) mode = ALL;
      else error->all(FLERR,"Illegal compute temp/em2 command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/em2 command");
  }

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempEM2::~ComputeTempEM2()
{
  delete [] id_bias;
  delete [] vector;

  delete [] conf_file;

  if (allocated) {
    memory->destroy(inertia);
    memory->destroy(quat_flag);
    memory->destroy(mem_flag);
    memory->destroy(prot_flag);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempEM2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(inertia,n+1,3,"fix:inertia");
  memory->create(quat_flag,n+1,"fix:quat_flag");
  memory->create(mem_flag,n+1,"fix:mem_flag");
  memory->create(prot_flag,n+1,"fix:protmem_flag");
}


/* ---------------------------------------------------------------------- */

void ComputeTempEM2::init()
{
  // error check

  avec = (AtomVecEM2 *) atom->style_match("em2");
  if (!avec)
    error->all(FLERR,"Compute temp/em2 requires atom style em2");

  // check that all particles are finite-size, no point particles allowed

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (tempbias) {
    int i = modify->find_compute(id_bias);
    if (i < 0) 
      error->all(FLERR,"Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
    if (tbias->tempflag == 0)
      error->all(FLERR,"Bias compute does not calculate temperature");
    if (tbias->tempbias == 0)
      error->all(FLERR,"Bias compute does not calculate a velocity bias");
    if (tbias->igroup != igroup)
      error->all(FLERR,"Bias compute group does not match compute group");
    tbias->init();
    tbias->setup();
    if (strcmp(tbias->style,"temp/region") == 0) tempbias = 2;
    else tempbias = 1;
  }

  if (!allocated) allocate();
  read_conf_file(conf_file);
}

/* ---------------------------------------------------------------------- */

void ComputeTempEM2::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  fix_dof = -1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempEM2::dof_compute()
{
  if (fix_dof) adjust_dof_fix();

  // 6 dof for 3d, 3 dof for 2d
  // which dof are included also depends on mode
  // assume full rotation of extended particles
  // user should correct this via compute_modify if needed

  double natoms = group->count(igroup);
  int nper;
  if (domain->dimension == 3) {
    if (mode == ALL) nper = 6;
    else nper = 3;
  } else {
    if (mode == ALL) nper = 3;
    else nper = 1;
  }
  dof = nper*natoms;

  // additional adjustments to dof

  if (tempbias == 1) {
    if (mode == ALL) dof -= tbias->dof_remove(-1) * natoms;

  } else if (tempbias == 2) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    tbias->dof_remove_pre();

    int count = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (tbias->dof_remove(i)) count++;
    int count_all;
    MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
    dof -= nper*count_all;
  }

  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempEM2::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double *inertia_i,*quat;
  double wbody[3];
  double rot[3][3];
  int itype;

  // sum translational and rotational energy for each particle
  // no point particles since divide by inertia

  double t = 0.0;

  if (mode == ALL) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        itype = type[i];

        if (rmass) {
          t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
        } else {
          t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * mass[itype];
        }

        quat = bonus[i].quat;
        inertia_i = inertia[itype];

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia_i[0];
        wbody[1] /= inertia_i[1];
        wbody[2] /= inertia_i[2];

        t += inertia_i[0]*wbody[0]*wbody[0] + inertia_i[1]*wbody[1]*wbody[1] + inertia_i[2]*wbody[2]*wbody[2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        itype = type[i];

        quat = bonus[i].quat;
        inertia_i = inertia[itype];

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia_i[0];
        wbody[1] /= inertia_i[1];
        wbody[2] /= inertia_i[2];

        t += inertia_i[0]*wbody[0]*wbody[0] + inertia_i[1]*wbody[1]*wbody[1] + inertia_i[2]*wbody[2]*wbody[2];
      }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic || tempbias == 2) dof_compute();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempEM2::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_vector != update->ntimestep) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  AtomVecEM2::Bonus *bonus = avec->bonus;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double *inertia_i,*quat;
  double wbody[3],t[6];
  double rot[3][3];
  double massone;
  int itype;

  // sum translational and rotational energy for each particle
  // no point particles since divide by inertia

  for (i = 0; i < 6; i++) t[i] = 0.0;

  if (mode == ALL) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        itype = type[i];

        if (rmass)
          massone = rmass[i];
        else
          massone = mass[itype];

        t[0] += massone * v[i][0]*v[i][0];
        t[1] += massone * v[i][1]*v[i][1];
        t[2] += massone * v[i][2]*v[i][2];
        t[3] += massone * v[i][0]*v[i][1];
        t[4] += massone * v[i][0]*v[i][2];
        t[5] += massone * v[i][1]*v[i][2];

        quat = bonus[i].quat;
        inertia_i = inertia[itype];

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia_i[0];
        wbody[1] /= inertia_i[1];
        wbody[2] /= inertia_i[2];

        // rotational kinetic energy

        t[0] += inertia_i[0]*wbody[0]*wbody[0];
        t[1] += inertia_i[1]*wbody[1]*wbody[1];
        t[2] += inertia_i[2]*wbody[2]*wbody[2];
        t[3] += inertia_i[0]*wbody[0]*wbody[1];
        t[4] += inertia_i[1]*wbody[0]*wbody[2];
        t[5] += inertia_i[2]*wbody[1]*wbody[2];
      }

  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        itype = type[i];

        quat = bonus[i].quat;
        inertia_i = inertia[itype];

        if (rmass)
          massone = rmass[i];
        else
          massone = mass[itype];

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia_i[0];
        wbody[1] /= inertia_i[1];
        wbody[2] /= inertia_i[2];

        // rotational kinetic energy

        t[0] += inertia_i[0]*wbody[0]*wbody[0];
        t[1] += inertia_i[1]*wbody[1]*wbody[1];
        t[2] += inertia_i[2]*wbody[2]*wbody[2];
        t[3] += inertia_i[0]*wbody[0]*wbody[1];
        t[4] += inertia_i[1]*wbody[0]*wbody[2];
        t[5] += inertia_i[2]*wbody[1]*wbody[2];
      }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempEM2::remove_bias(int i, double *v)
{
  if (tbias) tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempEM2::restore_bias(int i, double *v)
{
  if (tbias) tbias->restore_bias(i,v);
}

/* ---------------------------------------------------------------------- */

void ComputeTempEM2::read_conf_file(char *filename)
{
  int i;
  int ilo,ihi,count;
  int flag1, flag2, flag3;
  double I[3];
  int *type_flag;

  FILE *file;
  int narg=0;
  char ln[1024], *line, *arg[16];

  int n=atom->ntypes;
  memory->create(type_flag,n+1,"compute:type_flag");

  for (i=1;i<=n;i++) type_flag[i] = 0;

  file = fopen(filename,"r");
  if (!file) error->all(FLERR,"Compute temp/em2: Error opening configuration file!");

  while ( fgets ( ln, sizeof ln, file ) != NULL ) {
    line = trim(ln);

    if (line[0]=='#' || isEmptyString(line)) continue;

    // Split the parameter line into strings
    narg = 0;
    arg[narg] = strtok(line, " \t\n");
    while ( arg[narg]!=NULL ) {
      narg++;
      if (narg>=16) error->all(FLERR,"Compute temp/em2: Configuration file line is too long!");
      arg[narg] = strtok(NULL, " \t\n");
    }

    if (narg!=7) error->all(FLERR,"Compute temp/em2: Configuration file format error");

    I[0] = atof(arg[1]);
    I[1] = atof(arg[2]);
    I[2] = atof(arg[3]);
    flag1 = atoi(arg[4]);
    flag2 = atoi(arg[5]);
    flag3 = atoi(arg[6]);

    force->bounds(arg[0],atom->ntypes,ilo,ihi);

    if ((flag1!=0 && flag1!=1) || (flag2!=0 && flag2!=1) || (flag3!=0 && flag3!=1))
      error->all(FLERR,"Compute temp/em2: Configuration file format error");

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

    if (count==0) error->all(FLERR,"Compute temp/em2: Configuration file format error");
  }

  fclose(file);

  for (i=1;i<=n;i++) {
    if (!type_flag) error->all(FLERR,"Compute temp/em2: Parameters for all types must be set!");
  }

  memory->destroy(type_flag);
}

inline void ComputeTempEM2::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

char *ComputeTempEM2::ltrim(char *s)
{     
  while(isspace(*s)) s++;     
  return s; 
}  

char *ComputeTempEM2::rtrim(char *s)
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

char *ComputeTempEM2::trim(char *s)
{     
  return rtrim(ltrim(s));  
}

bool ComputeTempEM2::isEmptyString(char *str)
{
  int len = strlen(str);
  
  if (len==0) return true;
  
  for (int i=0;i<len;++i) {
    if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') return false;
  }

  return true;
}
