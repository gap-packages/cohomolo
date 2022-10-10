#include "defs.h"

#define NPT 32767
#define PSP 2000000
#define SVSP 200000
#define MP 1000
#define MB 80
#define MEXP 500

char mult, subgp, sgc, sgstr[2], inf0[80], inf1[80], inf2[80], inf3[80],
    inf4[80], outf[80], outft[80];
/* Defaults  inf1 gpname.spc
             inf2 gpname.dcr
             outf gpname.sc
         if -s set inf3 gpname.sg1 inf4 gpname.cr1.
         if -si set inf3 and inf4 are variable, and take values
         gpname.sgj and gpname.crj for 1<=j<=i.
*/
short perm[PSP], sv[SVSP], cp[2 * NPT], fpt[MEXP], orb[NPT + 1], intpow[MP],
    base[MB], lorb[MB], pno[MP], *pptr[MP], *svptr[MB], intno[MEXP],
    ngno[MEXP], power[MEXP], wt[MEXP], co[MEXP], rel[2 * MEXP],
    igno[2 * MEXP], pinv[NPT / 2], *sv2ptr[MB], d1[MEXP], d2[MEXP], pwt[MEXP],
    crep[NPT + 1], crepinv[NPT + 1], dcrep[NPT + 1], dcrepinv[NPT + 1],
    tp[2 * (NPT + 1)], mp = MP, mpt = NPT, mb = MB, mexp = MEXP;
int psp = PSP, svsp = SVSP;

int main(int argc, char * argv[])
{
  short arg;
  char  c, err, ondef;
  err = 0;
  ondef = 0;
  arg = 1;
  subgp = 0;
  mult = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    c = argv[arg][1];
    if (c == 's') {
      sgc = argv[arg][2];
      if (sgc == '\0') {
        subgp = 1;
        sgc = '1';
        ondef = 0;
      }
      else {
        subgp = sgc - '0';
        sgstr[1] = '\0';
        ondef = 1;
      }
      if (subgp <= 0 || subgp > 9) {
        err = 1;
        goto error;
      }
    }
    else if (c == 'm')
      mult = 1;
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
  strcpy(inf2, inf1);
  strcpy(outf, inf1);
  strcpy(outft, inf1);
  strcat(outft, "dcrt");
  if (subgp) {
    if (subgp == 1) {
      strcpy(inf3, inf1);
      strcpy(inf4, inf1);
      if (ondef == 0) {
        arg++;
        if (argc <= arg)
          strcat(inf3, "sg1");
        else
          strcat(inf3, argv[arg]);
        arg++;
        if (argc <= arg)
          strcat(inf4, "cr1");
        else
          strcat(inf4, argv[arg]);
      }
      else {
        strcat(inf4, "cr1");
        strcat(inf3, "sg1");
      }
    }
    else
      strcpy(inf0, inf1);
  }
  arg++;
  if (argc <= arg)
    strcat(outf, "sc");
  else
    strcat(outf, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(inf2, "dcr");
  else
    strcat(inf2, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(inf1, "spc");
  else
    strcat(inf1, argv[arg]);
  if (scprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr,
            "scrun [-s(n)] [-m] gpname (if -s [inf3(subgp)]) [inf4(cr)])\n");
    fprintf(stderr, "                [outf] [inf2(dcr)] [inf1].\n");
    exit(1);
  }
  exit(0);
}
