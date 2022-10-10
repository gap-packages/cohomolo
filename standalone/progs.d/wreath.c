#include "defs.h"

#define PSP 100000
#define MP 500
#define MB 80

char inf1[80], inf2[80], outf[80];
/* No defaults for filenames */
short perm[PSP], *pptr1[MP + 1], *pptr2[MP + 1], base[MB], mp = MP, mb = MB;
int   psp = PSP;
FILE *ip, *op;

int main(int argc, char * argv[])
{
  short npt1, npt2, i, j, k, l, m, n, np, np1, np2, *p, *q, mxp, *perm2, npto;
  char  err;
  int   quot;
  if (argc != 4) {
    err = 1;
    goto error;
  }
  else
    err = 0;
  strcpy(inf1, argv[1]);
  strcpy(inf2, argv[2]);
  strcpy(outf, argv[3]);
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    exit(1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt1, &np1, &m, &n);
  snl();
  if (m != 0)
    snl();
  if (n != 0)
    snl();
  quot = psp / npt1;
  if (quot > mp)
    quot = mp;
  mxp = quot;
  if (np1 > mxp) {
    fprintf(stderr, "np1 too big. Increase PSP and/or MP.\n");
    exit(1);
  }
  for (i = 1; i <= np1; i++)
    pptr1[i] = perm - 1 + (i - 1) * npt1;
  for (i = 1; i <= np1; i++) {
    p = pptr1[i];
    for (n = 1; n <= npt1; n++)
      fscanf(ip, "%hd", p + n);
    snl();
  }
  fclose(ip);
  psp -= (npt1 * np1);
  perm2 = perm + npt1 * np1;
  if ((ip = fopen(inf2, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf2);
    exit(1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt2, &np2, &m, &n);
  snl();
  if (m != 0)
    snl();
  if (n != 0)
    snl();
  mxp = psp / npt2;
  if (mxp > mp)
    mxp = mp;
  if (np2 > mxp) {
    fprintf(stderr, "np2 too big. Increase PSP and/or MP.\n");
    exit(1);
  }
  for (i = 1; i <= np2; i++)
    pptr2[i] = perm2 - 1 + (i - 1) * npt2;
  for (i = 1; i <= np2; i++) {
    p = pptr2[i];
    for (n = 1; n <= npt2; n++)
      fscanf(ip, "%hd", p + n);
    snl();
  }
  fclose(ip);
  /* Groups are now input  */
  op = fopen(outf, "w");
  npto = npt1 * npt2;
  fprintf(op, "%4d%4d%4d%4d\n", npto, np1 + np2, 0, 0);
  for (i = 1; i <= np1; i++) {
    p = pptr1[i];
    if (npto >= 1000) {
      for (j = 1; j <= npt1; j++)
        fprintf(op, "%5d", p[j]);
      for (j = npt1 + 1; j <= npto; j++)
        fprintf(op, "%5d", j);
    }
    else {
      for (j = 1; j <= npt1; j++)
        fprintf(op, "%4d", p[j]);
      for (j = npt1 + 1; j <= npto; j++)
        fprintf(op, "%4d", j);
    }
    fprintf(op, "\n");
  }
  for (i = 1; i <= np2; i++) {
    p = pptr2[i];
    if (npto >= 1000)
      for (j = 1; j <= npt2; j++)
        for (k = 1; k <= npt1; k++)
          fprintf(op, "%5d", (p[j] - 1) * npt1 + k);
    else
      for (j = 1; j <= npt2; j++)
        for (k = 1; k <= npt1; k++)
          fprintf(op, "%4d", (p[j] - 1) * npt1 + k);
    fprintf(op, "\n");
  }
  fclose(op);
error:;
  if (err) {
    printf("Usage: wreath inf1 inf2 outf.\n");
    exit(1);
  }
  exit(0);
}

int snl(void)
{
  while (getc(ip) != '\n')
    ;
}
