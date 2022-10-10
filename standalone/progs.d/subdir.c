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
  short npt1, npt2, i, l, m, n, np, np1, np2, *p, *q, mxp, *perm2, nb1, nb2,
      nb;
  char err;
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
  mxp = psp / npt1 - 1;
  if (mxp > mp)
    mxp = mp;
  if (np1 > mxp) {
    fprintf(stderr, "np1 too big. Increase PSP and/or MP.\n");
    exit(1);
  }
  for (i = 0; i <= np1; i++)
    pptr1[i] = perm - 1 + i * npt1;
  p = pptr1[0];
  for (i = 1; i <= npt1; i++)
    p[i] = i;
  for (i = 1; i <= np1; i++) {
    p = pptr1[i];
    for (n = 1; n <= npt1; n++)
      fscanf(ip, "%hd", p + n);
    snl();
  }
  fclose(ip);
  psp -= (npt1 * (np1 + 1));
  perm2 = perm + npt1 * (np1 + 1);
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
  mxp = psp / npt2 - 1;
  if (mxp > mp)
    mxp = mp;
  if (np2 > mxp) {
    fprintf(stderr, "np2 too big. Increase PSP and/or MP.\n");
    exit(1);
  }
  for (i = 0; i <= np2; i++)
    pptr2[i] = perm2 - 1 + i * npt2;
  p = pptr2[0];
  for (i = 1; i <= npt2; i++)
    p[i] = i;
  for (i = 1; i <= np2; i++) {
    p = pptr2[i];
    for (n = 1; n <= npt2; n++)
      fscanf(ip, "%hd", p + n);
    snl();
  }
  fclose(ip);
  /* Two groups are now input */
  op = fopen(outf, "w");
  printf("Do you wish to output any base points (y/n)?  ");
  if (getchar() == 'y') {
    printf("Input base points from group 1, preceded by their number.\n");
    scanf("%hd", &nb1);
    if (nb1 >= mb) {
      fprintf(stderr, "Too many base points. Increase MB.\n");
      exit(1);
    }
    for (i = 1; i <= nb1; i++)
      scanf("%hd", base + i);
    printf("Input base points from group 2, preceded by their number.\n");
    scanf("%hd", &nb2);
    nb = nb1 + nb2;
    if (nb >= mb) {
      fprintf(stderr, "Too many base points. Increase MB.\n");
      exit(1);
    }
    p = base + nb1 + 1;
    for (i = 1; i <= nb2; i++) {
      scanf("%hd", p);
      (*p) += npt1;
      p++;
    }
  }
  else
    nb = 0;
  printf("How many permutations do you wish to output?  ");
  scanf("%hd", &np);
  fprintf(op, "%4d%4d%4d%4d\n", npt1 + npt2, np, nb, 0);
  if (nb > 0)
    if (npt1 + npt2 >= 1000) {
      for (i = 1; i <= nb; i++)
        fprintf(op, "%5d", base[i]);
      fprintf(op, "\n");
    }
    else {
      for (i = 1; i <= nb; i++)
        fprintf(op, "%4d", base[i]);
      fprintf(op, "\n");
    }
  printf("Now input the %d pairs of perm nos. Perm 0 = identity.\n", np);
  for (i = 1; i <= np; i++) {
    scanf("%hd%hd", &m, &n);
    while (m > np1 || n > np2) {
      fprintf(stderr, "One of those nos is too big! Try again.\n");
      scanf("%hd%hd", &m, &n);
    }
    p = pptr1[m];
    q = pptr2[n];
    if (npt1 + npt2 >= 1000) {
      for (l = 1; l <= npt1; l++)
        fprintf(op, "%5d", p[l]);
      for (l = 1; l <= npt2; l++)
        fprintf(op, "%5d", q[l] + npt1);
    }
    else {
      for (l = 1; l <= npt1; l++)
        fprintf(op, "%4d", p[l]);
      for (l = 1; l <= npt2; l++)
        fprintf(op, "%4d", q[l] + npt1);
    }
    fprintf(op, "\n");
  }
error:
  if (err) {
    printf("Usage:    subdir inf1 inf2 outf.\n");
    exit(1);
  }
  exit(0);
}

int snl(void)
{
  while (getc(ip) != '\n')
    ;
}
