#include "defs.h"

#define MSP 50000
#define MV 500

short mv = MV, mspace[MSP], *vec[MV];
char  inf1[80], inf2[80], outf[80];
/* no defaults */
FILE *ip1, *ip2, *op;
int   msp = MSP;

int main(int argc, char * argv[])
{
  char  c, arg, err, eq, wedge;
  short i, j, k, l, n, prime, dim1, dim2, dimo, **m1, **m2, nmats, s;
  int   x;
  arg = 1;
  err = 0;
  eq = 0;
  wedge = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  if (argv[arg][0] == '-') {
    c = argv[arg][1];
    if (c == 'w') {
      wedge = 1;
      eq = 1;
    }
    else if (c == 'e')
      eq = 1;
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
  arg++;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcat(inf1, argv[arg]);
  if (eq == 0) {
    arg++;
    if (argc <= arg) {
      err = 1;
      goto error;
    }
    strcat(inf2, argv[arg]);
  }
  arg++;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcat(outf, argv[arg]);

  if ((ip1 = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s\n", inf1);
    exit(1);
  }
  if (eq == 0) {
    if ((ip2 = fopen(inf2, "r")) == 0) {
      fprintf(stderr, "Cannot open %s\n", inf1);
      exit(1);
    }
  }
  fscanf(ip1, "%hd%hd%hd", &prime, &dim1, &nmats);
  if (eq == 0) {
    fscanf(ip2, "%hd%hd%hd", &i, &dim2, &j);
    if (i != prime || j != nmats) {
      fprintf(stderr, "prime or nmats in two groups disagree.\n");
      exit(1);
    }
  }
  else
    dim2 = 0;
  if (dim1 + dim2 > mv) {
    fprintf(stderr, "Need more vectors. Increase MV.\n");
    exit(1);
  }
  if (dim1 * dim1 + dim2 * dim2 > msp) {
    fprintf(stderr, "Need more space. Increase MSP.\n");
    exit(1);
  }
  for (i = 0; i < dim1; i++)
    vec[i] = mspace + i * dim1 - 1;
  m1 = vec - 1;
  k = dim1 * dim1 - 1;
  if (eq == 0) {
    for (i = 0; i < dim2; i++)
      vec[i + dim1] = mspace + k + i * dim2;
    m2 = vec + dim1 - 1;
  }
  else {
    dim2 = dim1;
    m2 = m1;
  }
  dimo = wedge ? dim1 * (dim1 - 1) / 2 : dim1 * dim2;
  op = fopen(outf, "w");
  fprintf(op, "%3d%5d%3d\n", prime, dimo, nmats);
  for (n = 1; n <= nmats; n++) {
    for (i = 1; i <= dim1; i++)
      for (j = 1; j <= dim1; j++)
        fscanf(ip1, "%hd", m1[i] + j);
    if (eq == 0)
      for (i = 1; i <= dim2; i++)
        for (j = 1; j <= dim2; j++)
          fscanf(ip2, "%hd", m2[i] + j);
    if (wedge) {
      for (i = 1; i <= dim1; i++)
        for (j = 1; j < i; j++) {
          for (k = 1; k <= dim1; k++)
            for (l = 1; l < k; l++) {
              x = m1[i][k] * m1[j][l] - m1[i][l] * m1[j][k];
              s = x % prime;
              if (s < 0)
                s += prime;
              fprintf(op, "%4d", s);
            }
          fprintf(op, "\n");
        }
    }
    else {
      for (i = 1; i <= dim1; i++)
        for (j = 1; j <= dim2; j++) {
          for (k = 1; k <= dim1; k++)
            for (l = 1; l <= dim2; l++) {
              x = m1[i][k] * m2[j][l];
              s = x % prime;
              if (s < 0)
                s += prime;
              fprintf(op, "%4d", s);
            }
          fprintf(op, "\n");
        }
    }
  }
error:
  if (err) {
    printf("Usage: tens [-e] [-w] gpname inf1 (inf2) outf.\n");
    exit(1);
  }
  exit(0);
}
