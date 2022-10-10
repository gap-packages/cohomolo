#include "defs.h"

#define SPACE 2000000
#define CSPACE 1000000
#define WDL 20000
#define PTRSP 10000
#define CDPTRSP 10000
#define CPTRSP 500000
#define MB 80
#define MM 129
#define MPR 500
#define MARG 20000
/* This program uses output of extprun. Some of variables are the same */

char mult, inf0[80], inf1[80], inf2[80], outf[80], outft[80], inf3[80];
/* No defaults. inf1 is G inf2 H for corestriction H to G
   inf1 also takes values G.cp and G.rel.
   inf3 mats of cosetreps. (only if -m not set)
   outf always inf1.er. inf0 to remember gpname.
*/
short csp[CSPACE], ***coeff[MB], *cpsp[CPTRSP], **cdpsp[CDPTRSP];
short sp[SPACE], *psp[PTRSP], **imcos[MB], **mat[MM], **cpco[MB], lorb[MB],
    pinv[100], wdl = WDL, ptrsp = PTRSP, marg = MARG, cdptrsp = CDPTRSP,
               mb = MB, mpr = MPR;
int space = SPACE, cptrsp = CPTRSP, cspace = CSPACE;

int main(int argc, char * argv[])
{
  short arg;
  short err;
  mult = 0;
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  if (argv[arg][0] == '-') {
    if (argv[arg][1] != 'm') {
      err = 1;
      goto error;
    }
    mult = 1;
    arg++;
    if (argc <= arg) {
      err = 1;
      goto error;
    }
  }
  strcpy(inf0, argv[arg]);
  strcat(inf0, ".");
  strcpy(inf2, inf0);
  strcpy(inf3, inf0);
  arg++;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcat(inf0, argv[arg]);
  arg++;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcat(inf2, argv[arg]);
  if (mult == 0) {
    arg++;
    if (argc <= arg) {
      err = 1;
      goto error;
    }
    strcat(inf3, argv[arg]);
    strcat(inf3, "mat");
  }
  strcpy(outf, inf0);
  strcat(outf, ".er");
  strcpy(outft, inf0);
  strcat(outft, ".wl");
  if (crprog1() == -1)
    exit(1);
  if (crprog2() == -1)
    exit(1);
  ;
error:
  if (err) {
    fprintf(stderr,
            "Usage:    crrun [-m] gpname inf1 inf2 (if not -m inf3)\n");
    exit(1);
  }
  exit(0);
}
