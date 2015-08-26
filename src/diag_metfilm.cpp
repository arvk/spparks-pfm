/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_metfilm.h"
#include "app.h"
#include "app_metfilm.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{ZERO,FE,HYDROGEN,FEHMOL,VACANCY,SUP,SDOWN};      // same as DiagErbium
enum{ER,H,HE,VAC,EVENTS,ONE,TWO,THREE};

/* ---------------------------------------------------------------------- */

DiagMetfilm::DiagMetfilm(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"metfilm") != 0)
    error->all(FLERR,"Diag_style metfilm requires app_style metfilm");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Illegal diag_style metfilm command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style metfilm command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagMetfilm::~DiagMetfilm()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagMetfilm::init()
{
  appmetfilm = (AppMetfilm *) app;

  int none = appmetfilm->none;
  int ntwo = appmetfilm->ntwo;
  int nthree = appmetfilm->nthree;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"fe") == 0) which[i] = FE;
    else if (strcmp(list[i],"h") == 0) which[i] = H;
    else if (strcmp(list[i],"hmol") == 0) which[i] = FEHMOL;
    else if (strcmp(list[i],"vac") == 0) which[i] = VAC;
    else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all(FLERR,"Invalid value setting in diag_style metfilm");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all(FLERR,"Invalid value setting in diag_style metfilm");
      index[i] = n - 1;
    } else if (list[i][0] == 't') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all(FLERR,"Invalid value setting in diag_style metfilm");
      index[i] = n - 1;
    } else error->all(FLERR,"Invalid value setting in diag_style metfilm");
  }

  siteflag = 0;
  for (int i = 0; i < nlist; i++)
    if (which[i] == FE || which[i] == H || which[i] == FEHMOL || which[i] == VAC)
      siteflag = 1;

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagMetfilm::compute()
{
  int sites[5],ivalue;

  if (siteflag) {
    sites[FE] = sites[HYDROGEN] = sites[FEHMOL] = sites[VACANCY] = 0;
    int *element = appmetfilm->element;
    int nlocal = appmetfilm->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE) ivalue = sites[FE];
    else if (which[i] == H) ivalue = sites[HYDROGEN];
    else if (which[i] == FEHMOL) ivalue = sites[FEHMOL];
    else if (which[i] == VAC) ivalue = sites[VACANCY];
    else if (which[i] == EVENTS) ivalue = appmetfilm->nevents;
    else if (which[i] == ONE) ivalue = appmetfilm->scount[index[i]];
    else if (which[i] == TWO) ivalue = appmetfilm->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appmetfilm->tcount[index[i]];
    
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagMetfilm::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %d",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagMetfilm::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
