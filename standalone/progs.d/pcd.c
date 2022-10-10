#include "defs.h"

#define NPT 32767
#define PSP 2000000
#define SVSP 200000
#define MP 1000
#define MB 80
#define MEXP 500

/* MEXP-1 = maximal exponent */
char mult, conv, gens, inf1[80], inf2[80], outf1[80], outf2[80], outf3[80],
    outf4[80];
/* defaults inf1  gpname.sylp
            inf2  gpname.psg (no choice. Only if -c set)
            outf1 gpname.spc
            outf2 gpname.pcp
            outf3 gpname.psgwds (no choice. Only if -c set)
            outf4 gpname.pcpgens (no choice. Only -f -g set)
*/
short perm[PSP], sv[SVSP], cp[5 * NPT], fpt[MP], orb[NPT + 1], igno[MP],
    base[MB], lorb[MB], pno[MP], *pptr[MP], *svptr[MB], gno[MEXP], ngno[MEXP],
    power[MEXP], wt[MEXP], d1[MEXP], d2[MEXP], facord[MEXP], pinv[NPT / 2],
    rel[2 * MEXP], mp = MP, mpt = NPT, mb = MB, mexp = MEXP;
int svsp = SVSP, psp = PSP;

int main(int argc, char * argv[])
{
  short arg;
  char  err;
  mult = 0;
  conv = 0;
  gens = 0;
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    if (argv[arg][1] == 'm')
      mult = 1;
    else if (argv[arg][1] == 'c')
      conv = 1;
    else if (argv[arg][1] == 'g')
      gens = 1;
    else {
      err = 1;
      goto error;
    }
    arg++;
    if (argc <= arg) {
      err = 1;
      goto error;
    }
  }
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  strcpy(outf1, inf1);
  strcpy(outf2, inf1);
  strcpy(inf2, inf1);
  strcpy(outf3, inf1);
  strcpy(outf4, inf1);
  arg++;
  if (argc <= arg)
    strcat(inf1, "sylp");
  else
    strcat(inf1, argv[arg]);
  if (conv) {
    arg++;
    if (argc <= arg)
      strcat(inf2, "psg");
    else
      strcat(inf2, argv[arg]);
    arg++;
    if (argc <= arg)
      strcat(outf3, "psgwds");
    else
      strcat(outf3, argv[arg]);
  }
  if (gens)
    strcat(outf4, "pcpgens");
  arg++;
  if (argc <= arg)
    strcat(outf1, "spc");
  else
    strcat(outf1, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf2, "pcp");
  else
    strcat(outf2, argv[arg]);
  if (pcprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:    pcrun [-m] [-c] gpname\n");
    fprintf(stderr,
            "          [inf1] (if c [inf2] [outf3]) [outf1] [outf2].\n");
    exit(1);
  }
  exit(0);
}
