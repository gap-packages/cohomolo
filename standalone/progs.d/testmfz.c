#include "defs.h"

#define MSP 100000
#define MM 200
#define MV 1000
#define MDIM 100
char opt[5];
int  msp = MSP;
int  mv = MV, mm = MM, mdim = MDIM, dim, *spv, **spm, mspace[MSP], *vec[MV],
    **mat[MM], cp[500];
FILE *ip, *op;

int main(void)
{
  int  i, j, l, m, n, *p, **dp, maxv, maxm, ord;
  int  quot;
  char c, fault, flnm[80], f;
  printf("Do you wish to read matrices from a file(y/n)?  ");
  if ((c = getchar()) == 'y') {
    f = 1;
    seeknln();
  reenter:
    printf("Input filename:  ");
    scanf("%s", flnm);
    if ((ip = fopen(flnm, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", flnm);
      seeknln();
      goto reenter;
    }
    fscanf(ip, "%d%d", &dim, &n);
  }
  else {
    f = 0;
    printf("Input dimension:    ");
    scanf("%d", &dim);
  }
  if (dim > mdim) {
    fprintf(stderr, "dim too big.\n");
    exit(1);
  }
  seeknln();
  /* Set matrix and vector pointers as in mcd.c  */
  quot = msp / dim;
  if (quot >= mv)
    quot = mv - 1;
  maxv = quot;
  for (i = 1; i <= maxv; i++)
    vec[i] = mspace - 1 + (i - 1) * dim;
  maxm = (maxv - 10) / dim - 1;
  if (maxm >= mm)
    maxm = mm - 1;
  for (i = 0; i <= maxm; i++)
    mat[i] = vec + 10 + i * dim;
  dp = mat[0];
  p = *dp;
  for (i = 1; i <= dim * dim; i++)
    *(++p) = 0;
  for (i = 1; i <= dim; i++)
    dp[i][i] = 1;
  spm = mat[maxm];
  spv = spm[1];
  printf("10 free vectors and %d matrices can be defined.\n", maxm - 1);
  printf("Mat 0 is the identity.\n");
  if (f) {
    for (i = 1; i <= n; i++)
      readmat(mat[i]);
    fclose(ip);
  }
  while (1) {
    printf("Choose option. Print 'l' for list.\n");
    scanf("%s", opt);
    if (strcmp(opt, "l") == 0) {
      seeknln();
      printf("rv n  or  rm n        = read vec or mat no. n.\n");
      printf("pv  n  or  pm  n      = write vec or mat no. n.\n");
      printf("tr m n                = calc transpose m[n] of m[m].\n");
      printf("inv m n               = calc inverse m[n] of m[m].\n");
      printf("im l m n              = calc image v[m] of v[l] under m[n].\n");
      printf("c  l m n              = calc comm  v[m] of v[l] under m[n].\n");
      printf(
          "pr n l i(1)..i(l)     = calc prod m[n] of m[i(1)]...m[i(l)].\n");
      printf("sum n l i(1)...i(l)   = calc sum m[n] of m[i(1)]...m[i(l)].\n");
      printf("                        (mat nos may be negative for "
             "difference.)\n");
      printf("ch n i j k            = change m[n][i][j] to k.\n");
      printf("eq m n                = test equality of m[m],m[n].\n");
      printf("ord m n               = calc order of m[m] and put inverse in "
             "m[n].\n");
      printf("op filename           = output some matrices to filename.\n");
      printf("q                     = quit.\n");
    }
    else if (strcmp(opt, "rv") == 0) {
      scanf("%d", &n);
      seeknln();
      if (n <= 0 || n >= 10) {
        fprintf(stderr, "Invalid vector number %d\n", n);
        continue;
      }
      for (i = 1; i <= dim; i++)
        scanf("%d", vec[n] + i);
      seeknln();
    }
    else if (strcmp(opt, "rm") == 0) {
      scanf("%d", &n);
      seeknln();
      if (n <= 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      for (i = 1; i <= dim; i++)
        for (j = 1; j <= dim; j++)
          scanf("%d", mat[n][i] + j);
      seeknln();
    }
    else if (strcmp(opt, "pv") == 0) {
      scanf("%d", &n);
      seeknln();
      if (n <= 0 || n >= 10) {
        fprintf(stderr, "Invalid vector number %d\n", n);
        continue;
      }
      for (i = 1; i <= dim; i++)
        printf(" %3d", vec[n][i]);
      printf("\n");
    }
    else if (strcmp(opt, "pm") == 0) {
      scanf("%d", &n);
      seeknln();
      if (n < 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      for (i = 1; i <= dim; i++) {
        for (j = 1; j <= dim; j++)
          printf(" %3d", mat[n][i][j]);
        printf("\n");
      }
    }
    else if (strcmp(opt, "tr") == 0) {
      scanf("%d%d", &m, &n);
      seeknln();
      if (n <= 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      if (m < 0 || m >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", m);
        continue;
      }
      trans(mat[m], mat[n]);
    }
    else if (strcmp(opt, "inv") == 0) {
      scanf("%d%d", &m, &n);
      seeknln();
      if (n <= 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      if (m < 0 || m >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", m);
        continue;
      }
      inv(mat[m], mat[n]);
    }
    else if (strcmp(opt, "im") == 0 || strcmp(opt, "c") == 0) {
      scanf("%d%d%d", &l, &m, &n);
      seeknln();
      if (l <= 0 || l >= 10) {
        fprintf(stderr, "Invalid vector number %d\n", l);
        continue;
      }
      if (m <= 0 || m >= 10) {
        fprintf(stderr, "Invalid vector number %d\n", m);
        continue;
      }
      if (n < 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      if (strcmp(opt, "im") == 0)
        im(vec[l], vec[m], mat[n]);
      else
        comm(vec[l], vec[m], mat[n]);
    }
    else if (strcmp(opt, "pr") == 0) {
      scanf("%d%d", &n, cp);
      if (n <= 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        seeknln();
        continue;
      }
      fault = 0;
      for (i = 1; i <= *cp; i++) {
        scanf("%d", &m);
        if (m < 0 || m >= maxm) {
          fprintf(stderr, "Invalid matrix number %d\n", m);
          fault = 1;
          break;
        }
        cp[i] = m;
      }
      seeknln();
      if (fault)
        continue;
      prod(cp, mat[n]);
    }
    else if (strcmp(opt, "sum") == 0) {
      scanf("%d%d", &n, cp);
      if (n <= 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        seeknln();
        continue;
      }
      fault = 0;
      for (i = 1; i <= *cp; i++) {
        scanf("%d", &m);
        if (m < 0 || m >= maxm) {
          fprintf(stderr, "Invalid matrix number %d\n", m);
          fault = 1;
          break;
        }
        cp[i] = m;
      }
      seeknln();
      for (i = 1; i <= dim; i++)
        for (j = 1; j <= dim; j++) {
          m = 0;
          for (l = 1; l <= *cp; l++)
            if (cp[l] > 0)
              m += mat[cp[l]][i][j];
            else
              m -= mat[-cp[l]][i][j];
          mat[n][i][j] = m;
        }
    }
    else if (strcmp(opt, "ch") == 0) {
      scanf("%d%d%d%d", &n, &i, &j, &l);
      seeknln();
      if (n <= 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      if (i <= 0 || j <= 0 || i > dim || j > dim) {
        fprintf(stderr, "Invalid row or column no %d,%d\n", i, j);
        continue;
      }
      mat[n][i][j] = l;
    }
    else if (strcmp(opt, "eq") == 0) {
      scanf("%d%d", &m, &n);
      seeknln();
      if (n < 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      if (m < 0 || m >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", m);
        continue;
      }
      f = 1;
      for (i = 1; i <= dim; i++) {
        for (j = 1; j <= dim; j++)
          if (mat[m][i][j] != mat[n][i][j]) {
            f = 0;
            break;
          }
        if (f == 0)
          break;
      }
      if (f)
        printf("Matrices are equal.\n");
      else
        printf("Matrices differ in place  %d,%d.\n", i, j);
    }
    else if (strcmp(opt, "ord") == 0) {
      scanf("%d%d", &m, &n);
      seeknln();
      if (n <= 0 || n >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", n);
        continue;
      }
      if (m < 0 || m >= maxm) {
        fprintf(stderr, "Invalid matrix number %d\n", m);
        continue;
      }
      *cp = 1;
      cp[1] = m;
      prod(cp, mat[n]);
      prod(cp, spm);
      ord = 1;
      while (1) {
        f = 1;
        for (i = 1; i <= dim; i++) {
          for (j = 1; j <= dim; j++)
            if (mat[0][i][j] != spm[i][j]) {
              f = 0;
              break;
            }
          if (f == 0)
            break;
        }
        if (f) {
          printf("Order=%d\n", ord);
          break;
        }
        *cp = 1;
        cp[1] = maxm;
        prod(cp, mat[n]);
        ord++;
        *cp = 2;
        cp[1] = m;
        cp[2] = n;
        prod(cp, spm);
      }
    }
    else if (strcmp(opt, "op") == 0) {
      scanf("%s", flnm);
      op = fopen(flnm, "w");
      printf("Input matrix nos to be saved preceded by their number\n");
      scanf("%d", &n);
      fprintf(op, "%4d%4d\n", dim, n);
      fault = 0;
      for (i = 1; i <= n; i++) {
        scanf("%d", &j);
        if (j < 0 || j >= maxm) {
          fprintf(stderr, "Invalid matrix number %d\n", j);
          fault = 1;
          break;
        }
        printmat(mat[j]);
      }
      seeknln();
      if (fault) {
        fclose(op);
        unlink(flnm);
        continue;
      }
      fclose(op);
    }
    else if (strcmp(opt, "q") == 0) {
      seeknln();
      break;
    }
    else {
      seeknln();
      printf("Invalid option.\n");
    }
  }
  exit(0);
}

void seeknln(void)
{
  while (getchar() != '\n')
    ;
}
