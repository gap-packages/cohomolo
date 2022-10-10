#include "defs.h"

char inf[80], outf[80];
/* Defaults inf=gpname.inperm, outf=gpname.inmat */
FILE *ip, *op;

int main(int argc, char * argv[])
{
  short npt, np, i, j, n, im, im1, im2, start, arg, prime, dim, const1, opn;
  char  upper, lower, heart, err, c;
  arg = 1;
  err = 0;
  upper = 0;
  lower = 0;
  heart = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  if (argv[arg][0] == '-') {
    c = argv[arg][1];
    if (c == 'u')
      upper = 1;
    else if (c == 'l')
      lower = 1;
    else if (c == 'h')
      heart = 1;
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
  strcpy(inf, argv[arg]);
  strcat(inf, ".");
  strcpy(outf, inf);
  arg++;
  if (argc <= arg)
    strcat(inf, "inperm");
  else
    strcat(inf, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf, "inmat");
  else
    strcat(outf, argv[arg]);
  if ((ip = fopen(inf, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf);
    exit(1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &np, &i, &j);
  seeknln();
  if (i > 0)
    seeknln();
  if (j != 0)
    seeknln();
  printf("Input prime!    ");
  scanf("%hd", &prime);
  if (heart && npt % prime != 0) {
    fprintf(stderr,
            "Cannot compute heart when prime does not divide degree.\n");
    exit(1);
  }
  if (upper || lower) {
    dim = npt - 1;
    start = 2;
  }
  else if (heart) {
    dim = npt - 2;
    start = 3;
  }
  else {
    dim = npt;
    start = 1;
  }
  op = fopen(outf, "w");
  fprintf(op, "%4d%4d%4d\n", prime, dim, np);

  for (n = 1; n <= np; n++) {
    if (lower || upper || heart)
      fscanf(ip, "%hd", &im1);
    if (upper)
      im1 = 0;
    if (heart)
      fscanf(ip, "%hd", &im2);
    for (i = start; i <= npt; i++) {
      fscanf(ip, "%hd", &im);
      if ((upper && im == 1) || (heart && im == 2))
        const1 = prime - 1;
      else if (heart && im1 == 2)
        const1 = 1;
      else
        const1 = 0;
      for (j = start; j <= npt; j++) {
        opn = (j == im)    ? (const1 + 1) % prime
              : (j == im1) ? (const1 + prime - 1) % prime
                           : const1;
        if (prime < 10)
          fprintf(op, "%2d", opn);
        else if (prime < 100)
          fprintf(op, "%3d", opn);
        else
          fprintf(op, "%4d", opn);
      }
      fprintf(op, "\n");
    }
    seeknln();
    fprintf(op, "\n");
  }
error:
  if (err) {
    fprintf(stderr,
            "Usage:    permmat [-u] [-l] [-h] gpname [inf] [outf].\n");
    exit(1);
  }
  exit(0);
}

void seeknln(void)
{
  while (getc(ip) != '\n')
    ;
}
