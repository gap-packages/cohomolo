#include "defs.h"

#define NPT 32767
#define PSP 1000000
#define SVSP 400000
#define MP 1000
#define MB 80
#define SPACE 500000
#define MEXP 500
#define RISP 32767
#define PAR1 5
#define PAR2 5
#define PAR3 5
#define PAR4 500
/* SPACE is for cosetrep perms and various arrays defined by pointers
   RISP is size of one of these
   MEXP = max no of perms that can be stored as coset reps,
   PAR1,..,PAR4 are parameters used in search for Sylow-group.
*/

char cent, sym, opt, hgst, nop, chpar, syl, nonb[NPT], inf1[80], inf2[80],
    inf3[80], outf1[80], outf2[80];
/* If -s not called, identical to normd.c.
   If -s called inf1 defaults to gpname.sg
           and  inf2 to gpname.sylp
         outf1 becomes gpname.temp, used to store subgroup of sylp-group
         found so far,which is then used as inf2 input to normalizer program.
         This is removed when no longer needed.
*/
short mp = MP, mexp = MEXP, mb = MB - 1, mnpt = NPT, risp = RISP, par1 = PAR1,
      par2 = PAR2, par3 = PAR3, par4 = PAR4, perm[PSP], sv[SVSP],
      cp[10 * NPT], orb[1 + NPT], gbase[NPT], hbase[MB], obase[MB], nbase[MB],
      lorbg[NPT], lorbn[NPT], lorbh[MB], ntno[NPT], reg[NPT + 1], ntorno[NPT],
      tsv1[NPT + 1], tsv2[NPT + 1], tsv3[NPT + 1], genorb[NPT + 1],
      expcp[NPT], fp[MP], pno[MP / 2], start[NPT + 1], space[SPACE], ipno[MP],
      endorno[NPT], *pptr[MP], *svgptr[MB], *svhptr[MB], *svnptr[NPT],
      *intorb[MB], *horno[MB], *hlorb[MB], *expptr[MEXP], *imorno[MB],
      *imlorb[MB], *orbperm[MB], *deftime[MB], *regsv[MB], orep[NPT + 1],
      prime;
int psp = PSP, sp = SPACE, svsp = SVSP;

int main(int argc, char * argv[])
{
  short arg, x;
  char  c, err;
  err = 0;
  arg = 1;
  hgst = 0;
  opt = 0;
  cent = 0;
  nop = 0;
  chpar = 0;
  syl = 0;
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
    else if (c == 's')
      syl = 1;
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
  if (syl) {
    hgst = 0;
    nop = 0;
    if (cent)
      cent = 0;
    if (opt)
      chpar = 1;
    strcpy(inf1, argv[arg]);
    strcat(inf1, ".");
    strcpy(inf2, inf1);
    strcpy(outf1, inf1);
    strcat(outf1, "temp");
    arg++;
    if (argc <= arg)
      strcat(inf1, "sg");
    else
      strcat(inf1, argv[arg]);
    arg++;
    if (argc <= arg)
      strcat(inf2, "sylp");
    else
      strcat(inf2, argv[arg]);
    strcpy(inf3, inf1);
    x = 0;
    while ((x = sylprog(-x)) > 0) {
      strcpy(inf1, outf1);
      if (nprg1() == -1)
        exit(1);
      if (nprg2() == -1)
        exit(1);
      while (sylprog(x) > 0) {
        if (nprg1() == -1)
          exit(1);
        if (nprg2() == -1)
          exit(1);
      }
      if (x == -1)
        exit(1);
      strcpy(inf1, inf3);
      unlink(outf1);
    }
    if (x == -1)
      exit(1);
  }
  else {
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
    if (nprg1() == -1)
      exit(1);
    if (nprg2() == -1)
      exit(1);
  }
error:
  if (err) {
    fprintf(stderr, "Usage:    sylnorm [-s] [-h] [-n] [-c] [-o] gpname\n");
    fprintf(stderr, "          inf1 inf2 (inf3) [outf1] [outf2].\n");
    exit(1);
  }
  exit(0);
}
