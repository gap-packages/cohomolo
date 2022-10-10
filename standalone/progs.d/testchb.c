#include "defs.h"
#include "permfns.h"

#define NPT 4000
#define MB 80
#define MP 800
#define PSP 200000
#define SVSP 50000

char  gpname[80];
short npt, nf, perm[PSP], sv[SVSP], cp[5 * NPT], orb[NPT + 1], base[MB],
    lorb[MB], pno[MP / 2], fp[MP], tsv1[NPT + 1], tsv2[NPT + 1],
    tsv3[NPT + 1], orep[NPT + 1], *pptr[MP], *svptr[MB], mp = MP, mb = MB - 1,
                                                         mnpt = NPT;
int   psp = PSP, svsp = SVSP;
FILE *ip, *op;

int main(int argc, char * argv[])
{
  short i, j, k, nperms, nb, st, pt, pos, mxp, mnb, arg;
  char  err, inf[80], outf[80];
  int   quot;
  /* Defaults: inf  gpname.sg
              outf  gpname.newbas
  */
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
    strcat(inf, "sg");
  else
    strcat(inf, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf, "newbas");
  else
    strcat(outf, argv[arg]);
  if ((ip = fopen(inf, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf);
    exit(1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &k);
  if (npt > mnpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    exit(1);
  }
  if (k <= 2) {
    fprintf(stderr, "Wrong input format.\n");
    exit(1);
  }
  quot = psp / (npt + 1);
  if (quot > mp)
    quot = mp;
  mxp = quot;
  if (mxp % 2 == 1)
    mxp--;
  quot = svsp / npt;
  if (quot > mb)
    quot = mb;
  mnb = quot;
  if (mnb >= npt)
    mnb = npt - 1;
  mb = mnb;
  if (nb > mb) {
    fprintf(stderr, "nb too big. Increase SVSP (or MB).\n");
    exit(1);
  }
  if (2 * nperms + 1 >= mxp) {
    fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
    exit(1);
  }
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + i * (npt + 1) - 1;
  for (i = 1; i <= mnb; i++)
    svptr[i] = sv + (i - 1) * npt - 1;
  readbaselo(nb, base, lorb);
  readpsv(0, nb, nperms, svptr);
  fclose(ip);
  for (i = 0; i < mxp; i += 2)
    fp[i] = i + 2;
  fp[mxp - 2] = -1;
  fp[2 * nperms - 2] = -1;
  nf = 2 * nperms;
  st = 0;
  while (1) {
    printf("Base and lorb are:\n");
    for (i = 1; i <= nb; i++)
      printf("%4d", base[i]);
    printf("\n");
    for (i = 1; i <= nb; i++)
      printf("%4d", lorb[i]);
    printf("\n");
    printf("Input new base and position (0 0 to quit):    ");
    scanf("%hd%hd", &pt, &pos);
    if (pt == 0)
      break;
    if ((i = intbase(pt, pos, &st, &nb, base, lorb, svptr)) == -1)
      exit(1);
  }
  *pno = 0;
  for (i = st; i != -1; i = fp[i]) {
    (*pno)++;
    pno[*pno] = i;
  }
  op = fopen(outf, "w");
  fprintf(op, "%4d%4d%4d%4d\n", npt, *pno, nb, 3);
  printbaselo(nb, base, lorb);
  printpsv(nb, pno, svptr);
error:
  if (err) {
    fprintf(stderr, "Usage:    chbrun gpname [inf] [outf].\n");
    exit(1);
  }
  exit(0);
}
