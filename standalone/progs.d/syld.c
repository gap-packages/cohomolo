#include "defs.h"

#define NPT 32767
#define PSP 1000000
#define SVSP 100000
#define MP 1000
#define MB 80
#define EXPSP 100000
#define MEXP 5000
/* EXPSP = space for stored cosetrep perms, MEXP = max no of these */

char inf[80], outf[80];
/* Defaults: inf gpname.outperm
            outf gpname.sylp
*/
short mp = MP, mxexp = MEXP, mb = MB - 1, mnpt = NPT, perm[PSP], svp2[SVSP],
      svg2[SVSP], cp[5 * NPT], orb[NPT + 1], expperm[EXPSP], adpt[MB],
      expcp[MB], lorbg[MB], lorbp[MB], lorbdef[MB], base[MB], ntfpt[MB],
      ntbpt[MB], start[MB + 1], invbase[NPT + 1], pno[MP / 2], facord[MP],
      orno[NPT + 1], lporb[NPT + 1], deftime[NPT + 1], orbperm[NPT + 1],
      genorb[NPT + 1], *pptr[MP], *svgptr[MB], *svpptr[MB], *expptr[MEXP];
int psp = PSP, expsp = EXPSP, svsp = SVSP;

int main(int argc, char * argv[])
{
  short arg;
  char  err;
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcpy(inf, argv[arg]);
  strcat(inf, ".");
  strcpy(outf, inf);
  arg++;
  if (argc <= arg)
    strcat(inf, "outperm");
  else
    strcat(inf, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf, "sylp");
  else
    strcat(outf, argv[arg]);
  if (sylprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:    sylrun gpname [inf] [outf].\n");
    exit(1);
  }
  exit(0);
}
