#include "defs.h"

#define NPT 32767
#define PSP 1000000
#define SVSP 100000
#define MP 500
#define MB 80
#define SPACE 1000000
#define MEXP 4000
#define RISP 32767
/* SPACE is for cosetrep perms and various arrays defined by pointers
   MEXP is max no of stored cosetrep perms (as in syld.c)
   RISP is size of space set aside for reular orbit info
       (this cannot be set bigger than 2^15-1).
*/

char cent, sym, opt, hgst, nop, nonb[NPT], inf1[80], inf2[80], inf3[80],
    outf1[80], outf2[80];
/*   inf1 =gpname.?1 (can be sym. No default)
     inf2 =gpname.?2 (must be subgroup of ?1. No default)
     inf3 =gpname.?3 (Usually not used. Small generating set for ?2.
                      No default)
     Defaults:  outf1 gpname.norm (or cent, if -c set)
                outf2 gpname.ng (only if -n called)
*/
short mp = MP, mexp = MEXP, mb = MB - 1, mnpt = NPT, risp = RISP, prime = 0,
      perm[PSP], sv[SVSP], cp[10 * NPT], orb[1 + NPT], gbase[NPT], hbase[MB],
      obase[MB], nbase[MB], lorbg[NPT], lorbn[NPT], lorbh[MB], ntno[NPT],
      reg[NPT], ntorno[NPT], tsv1[NPT + 1], tsv2[NPT + 1], tsv3[NPT + 1],
      genorb[NPT + 1], expcp[NPT], fp[MP], pno[MP / 2], start[NPT + 1],
      space[SPACE], ipno[MP], endorno[NPT], *pptr[MP], *svgptr[MB],
      *svhptr[MB], *svnptr[NPT], *intorb[MB], *horno[MB], *hlorb[MB],
      *expptr[MEXP], *imorno[MB], *imlorb[MB], *orbperm[MB], *deftime[MB],
      *regsv[MB], orep[NPT + 1];
int psp = PSP, sp = SPACE, svsp = SVSP;

int main(int argc, char * argv[])
{
  short arg, r;
  char  c, err;
  err = 0;
  arg = 1;
  hgst = 0;
  opt = 0;
  cent = 0;
  nop = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    c = argv[arg][1];
    if (c == 'h')
      hgst = 1;
    else if (c == 'o')
      opt = 1;
    else if (c == 'c')
      cent = 1;
    else if (c == 'n')
      nop = 1;
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
  if (argc <= arg + 2) {
    err = 1;
    goto error;
  }
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  strcpy(inf2, inf1);
  strcpy(inf3, inf1);
  strcpy(outf1, inf1);
  strcpy(outf2, inf1);
  if (strcmp(argv[arg + 1], "sym") == 0)
    sym = 1;
  else {
    strcat(inf1, argv[arg + 1]);
    sym = 0;
  }
  strcat(inf2, argv[arg + 2]);
  arg += 3;
  if (hgst) {
    if (argc <= arg) {
      err = 1;
      goto error;
    }
    strcat(inf3, argv[arg]);
    arg++;
  }
  if (argc <= arg) {
    if (cent)
      strcat(outf1, "cent");
    else
      strcat(outf1, "norm");
  }
  else {
    strcat(outf1, argv[arg]);
    arg++;
  }
  if (nop) {
    if (argc <= arg)
      strcat(outf2, "ng");
    else
      strcat(outf2, argv[arg]);
  }
  if (sym)
    printf("G is the symmetric group.\n");
  if (cent)
    printf("Calculation of C(H) ^ G.\n");
  else
    printf("Calculation of N(H) ^ G.\n");
  if ((r = nprg1()) == -1)
    exit(1);
  if (r == 1)
    exit(2); /* This means H = G */
  if ((r = nprg2()) == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:    normrun [-h] [-n] [-c] [-o] ");
    fprintf(stderr, "gpname inf1 inf2 (inf3) [outf1] [outf2].\n");
    exit(1);
  }
  exit(0);
}
