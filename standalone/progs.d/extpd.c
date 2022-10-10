#include "defs.h"

#define PSP 500000
#define SPACE 500000
#define PTRSP 10000
#define CSPACE 2000000
#define CPTRSP 500000
#define CDPTRSP 10000
#define SVSP 200000
#define NPT 32767
#define MP 2000
#define MB 80
#define MWDL 5000
#define MWL2 1000
#define MLWDL 32767
#define MSP 32767
#define MV 1000
#define MM 200
#define MDIM 101
#define MPR 2000

char     mult, inf1[80], inf2[80], inf3[80], outf[80], outft[80];
short ***scoeff[MB], **cdpsp[CDPTRSP], *cpsp[CPTRSP], csp[CSPACE],
    rwd[MLWDL + MDIM], wd1[MWDL + MDIM], wd2[MWDL + MDIM], wd3[MWDL + MDIM],
    wd4[MWDL + MDIM];
int   cspace = CSPACE, psp = PSP, space = SPACE, cptrsp = CPTRSP, svsp = SVSP;
short perm[PSP], sv[SVSP], orb[NPT + 1], imsp[SPACE], *ptsp[PTRSP],
    **simcos[MB], base[MB], lorb[MB], pno[MP / 2], *pptr[MP], *svptr[MB],
    cord[MDIM], pinv[MPR], invg[MP], rno[MB], cp[MLWDL + MDIM], mspace[MSP],
    *vec[MV], **mat[MM], gno[MB],
    mp = MP, cdptrsp = CDPTRSP, ptrsp = PTRSP, mm = MM, msp = MSP, mv = MV,
    mwdl = MWDL, mlwdl = MLWDL, mpr = MPR, mpt = NPT, mb = MB - 1,
    mdim = MDIM - 1, mwl2 = MWL2;

int main(int argc, char * argv[])
{
  short arg;
  char  err;
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
      arg = 1;
      goto error;
    }
  }
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  arg++;
  if (argc <= arg)
    strcat(inf1, "sg");
  else
    strcat(inf1, argv[arg]);
  strcpy(inf2, inf1);
  strcat(inf2, ".er");
  strcpy(inf3, inf1);
  strcat(inf3, "mat");
  strcpy(outf, inf1);
  strcat(outf, ".ep");
  strcpy(outft, inf1);
  strcat(outft, ".gc");
  if (extpprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:   extprun [-m] gpname [inf1]\n");
    exit(1);
  }
  exit(0);
}
