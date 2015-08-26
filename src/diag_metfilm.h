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

#ifdef DIAG_CLASS
DiagStyle(metfilm,DiagMetfilm)

#else

#ifndef SPK_DIAG_METFILM_H
#define SPK_DIAG_METFILM_H

#include "diag.h"

namespace SPPARKS_NS {

class DiagMetfilm : public Diag {
 public:
  DiagMetfilm(class SPPARKS *, int, char **);
  ~DiagMetfilm();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppMetfilm *appmetfilm;
  int nlist;
  char **list;
  int *which,*index,*ivector;
  int siteflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag_style metfilm requires app_style metfilm

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Invalid value setting in diag_style metfilm

UNDOCUMENTED

*/
