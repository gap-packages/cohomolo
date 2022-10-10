#include "defs.h"


char inf1[80], inf2[80], outf[80];
/* No defaults */
FILE *ip1, *ip2, *op;

int main(int argc, char * argv[])
{
  short prime, dim1, dim2, nmats, i, j, k, n, err;
  if (argc <= 4) {
    err = 1;
    goto error;
  }
  err = 0;
  strcpy(inf1, argv[1]);
  strcat(inf1, ".");
  strcpy(inf2, inf1);
  strcpy(outf, inf1);
  strcat(inf1, argv[2]);
  strcat(inf2, argv[3]);
  strcat(outf, argv[4]);
  if ((ip1 = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    exit(1);
  }
  if ((ip2 = fopen(inf2, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf2);
    exit(1);
  }
  fscanf(ip1, "%hd%hd%hd", &prime, &dim1, &nmats);
  fscanf(ip2, "%hd%hd%hd", &i, &dim2, &k);
  if (prime != i || nmats != k) {
    fprintf(stderr, "Two groups have different primes or nmats.\n");
    exit(1);
  }
  op = fopen(outf, "w");
  fprintf(op, "%3d%5d%3d\n", prime, dim1 + dim2, nmats);

  for (i = 1; i <= nmats; i++) {
    for (j = 1; j <= dim1; j++) {
      for (k = 1; k <= dim1; k++) {
        fscanf(ip1, "%hd", &n);
        fprintf(op, "%4d", n);
      }
      for (k = 1; k <= dim2; k++)
        fprintf(op, "%4d", 0);
      fprintf(op, "\n");
    }
    for (j = 1; j <= dim2; j++) {
      for (k = 1; k <= dim1; k++)
        fprintf(op, "%4d", 0);
      for (k = 1; k <= dim2; k++) {
        fscanf(ip2, "%hd", &n);
        fprintf(op, "%4d", n);
      }
      fprintf(op, "\n");
    }
  }

error:
  if (err) {
    printf("Usage:    sum gpname inf1 inf2 outf.\n");
    exit(1);
  }
  exit(0);
}
