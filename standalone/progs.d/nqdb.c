#include "defs.h"

#define MNG 500001
#define MEXP 500001
#define RSP 100000000
#define PTRSP 10000000
#define MCL 31
#define MSP 500000
#define MV 25000
#define MM 60
#define MWL 5000
#define MARG 500000
#define MPR 500

/* MEXP-1=max exponent, MNG-1=max no of new gens introduced, RSP=basic space,
   PTRSP=basic pointer space, MCL=max class, MSP,MV,MM as in mcd.c,
   MARG=safety margin for collections (effect is unpredictable if this is
   too small - I think it should be at least 2*MEXP!), MPR=max prime, MWL=max
   length of words for matrix computations (done in array cp) (no warnings
   given if this is too small!).
 */

char inf0[80], inf1[80], inf2[80], inf3[80], inf4[80], inf[80], outf0[80],
    outf1[80], outf2[80], outfd[80], outcopy[80], act, ch1, crel, cfm, gap;
int mv = MV, mm = MM, mexp = MEXP, mng = MNG, mcl = MCL, facexp, tails, depth,
    no, prime, exp, nng, class, *rpf, *rpb, *eexpnt, *enexpnt, **pcb, intexp,
    stage, dim, onng, **opcb, **npcb, *nd1, *nd2, *spv, **spm, rel[RSP],
    expnt[MEXP], nexpnt[MNG], prvec[MNG], pinv[MPR], wt[MEXP + MNG],
    d1[MEXP + MNG], d2[MEXP + MNG], *pcptr[PTRSP], **powptr[MEXP],
    **comptr[MEXP], *sspc[MCL], *sspf[MCL], sgen[MCL], sex[MCL], spgen[MCL],
    spex[MCL], spugen[MCL], dpth[MEXP + MNG], sd1[MEXP], sd2[MEXP], swt[MEXP],
    mspace[MSP], *vec[MV], **mat[MM], cp[MWL];
int rsp = RSP, msp = MSP, ptrsp = PTRSP, wsp, marg = MARG;

/*  act=1 if -a set, ch1=1 if -1 set, crel=1 if -c set,
    inf0 (if act) = pcp output from previous run of nqrun,
    inf1 = current pcp input file (variable),
    inf2 = matrices of min gen set for P,
    inf3 (if act) = output from scrun,
    inf4 (if act) = matrices of dcreps,
    inf  remembers group name,
    outf0 = pcp output file used in comp of H^2(P,M) or H^2(Q,M),
    outf1 = current pcp output file (variable),
    outf2 = matrices of pcp gens of P,
    outfd = output of dimension of cohomolgy group for GAP
*/

int main(int argc, char * argv[])
{
  int  arg, n;
  char err;
  act = 0;
  ch1 = 0;
  cfm = 0;
  err = 0;
  arg = 1;
  crel = 0;
  gap = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    if (argv[arg][1] == 'a')
      act = 1;
    else if (argv[arg][1] == '1')
      ch1 = 1;
    else if (argv[arg][1] == 'c')
      crel = 1;
    else if (argv[arg][1] == 'f')
      cfm = 1;
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
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  strcpy(outf1, inf1);
  strcpy(inf2, inf1);
  strcpy(outf2, inf1);
  strcpy(outf0, inf1);
  strcpy(inf, inf1);
  strcpy(outfd, inf1);
  strcat(outfd, "cdim");
  if (ch1) {
    strcpy(outcopy, inf1);
    strcat(outcopy, "copy");
  }
  if (act) {
    strcpy(inf3, inf1);
    strcpy(inf4, inf1);
    strcpy(inf0, inf1);
    arg++;
    if (argc <= arg)
      strcat(inf3, "sc");
    else
      strcat(inf3, argv[arg]);
    arg++;
    if (argc <= arg)
      strcat(inf4, "dcrmat");
    else
      strcat(inf4, argv[arg]);
  }
  arg++;
  if (argc <= arg) {
    if (act) {
      if (ch1)
        strcat(inf0, "ch1");
      else
        strcat(inf0, "ch2");
    }
    else
      strcat(inf1, "pcp");
  }
  else {
    if (act)
      strcat(inf0, argv[arg]);
    else
      strcat(inf1, argv[arg]);
  }
  arg++;
  if (argc <= arg) {
    if (act)
      strcat(outf0, "octemp");
    else if (ch1)
      strcat(outf0, "ch1");
    else
      strcat(outf0, "ch2");
  }
  else
    strcat(outf0, argv[arg]);
  if (act == 0) {
    arg++;
    if (argc <= arg)
      strcat(inf2, "pgmat");
    else
      strcat(inf2, argv[arg]);
  }
  arg++;
  if (argc <= arg)
    strcat(outf2, "pmats");
  else
    strcat(outf2, argv[arg]);
  n = nqprog();
  if (n == -1)
    exit(1);
  if (n == 2)
    exit(2);
  /* This means cohomology group has become trivial */
  exit(0);
error:
  if (err) {
    fprintf(stderr, "Usage:    nqrun");
    fprintf(stderr, " [-1] [-c] [-a] [-f] gpname (if a [inf3] [inf4] )\n");
    fprintf(stderr, "      [inf1] [outf1] (if not a [inf2] [outf2]).\n");
    exit(1);
  }
}
