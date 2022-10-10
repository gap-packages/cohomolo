#include "defs.h"

#define MNG 1001
#define MEXP 41
#define RSP 100000
#define PTRSP 10000
#define WKSP 2500
/* MNG-1=max no of gens in multiplier, MEXP-1=max no of gens of P,
   WKSP=safety margin for various words
   RSP=amount of space for PCP relations
   PTRSP= amount of space for pointers to PCP relations.
*/

char gap, act, crel, ims, inf0[80], inf1[80], inf2[80], outf[80], outfm[80];
/* Defaults: inf1=gpname.pcp if not -a
                 =gpname.cov if -a
             outf=gpname.cov
             outfm=gpname.mult
             inf2=gpname.sc
             inf0 used to remember gpname.
*/
short mng = MNG - 1, mexp = MEXP - 1, wksp = WKSP, prime, exp, nng, class,
      *rpf, *rpb, *eexpnt, *enexpnt, **pcb, mnng, mord, intexp, rel[RSP],
      expnt[MEXP], nexpnt[MNG], cord[MNG], wt[MEXP], d1[MEXP + MNG],
      d2[MEXP + MNG], *pcptr[PTRSP], **powptr[MEXP], **comptr[MEXP],
      *sspc[MEXP], *sspf[MEXP], sgen[MEXP], sex[MEXP], spgen[MEXP],
      spex[MEXP], spugen[MEXP], *intg[MEXP], *imintg[MEXP], *tlintg[MEXP];
int rsp = RSP, ptrsp = PTRSP;

int main(int argc, char * argv[])
{
  short arg;
  char  err;
  ims = 0;
  act = 0;
  crel = 0;
  gap = 0;
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    if (argv[arg][1] == 'a')
      act = 1;
    else if (argv[arg][1] == 'c')
      crel = 1;
    else if (argv[arg][1] == 'i')
      ims = 1;
    else if (argv[arg][1] == 'g')
      gap = 1;
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
  strcpy(inf0, argv[arg]);
  strcat(inf0, ".");
  strcpy(outf, inf0);
  strcpy(outfm, inf0);
  strcat(outfm, "mult");
  strcpy(inf1, inf0);
  if (act) {
    strcpy(inf2, inf1);
    arg++;
    if (argc <= arg)
      strcat(inf2, "sc");
    else
      strcat(inf2, argv[arg]);
  }
  arg++;
  if (argc <= arg)
    if (act || crel)
      strcat(inf1, "cov");
    else
      strcat(inf1, "pcp");
  else
    strcat(inf1, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf, "cov");
  else
    strcat(outf, argv[arg]);
  if (nqmprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:  nqmrun");
    fprintf(stderr, " [-a] [-c] [-i] [-g] gpname ([inf2]) [inf1] [outf].\n");
    exit(1);
  }
  exit(0);
}
