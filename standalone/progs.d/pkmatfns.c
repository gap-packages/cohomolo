#include "defs.h"

extern char  prime, **mat[], cvec[], pinv[], aut;
extern short dim, svec[], maxnull;
extern FILE *ip, *op;

int trans(char ** a, char ** b)
/* The transpose of matrix a is written into matrix b.
   Externals: dim.
*/
{
  short i, j;
  for (i = 1; i <= dim; i++) {
    b[i][i] = a[i][i];
    for (j = 1; j < i; j++) {
      b[i][j] = a[j][i];
      b[j][i] = a[i][j];
    }
  }
}

int copy(char ** a, char ** b)
/* copies matrix a to matrix b.  Externals:  dim  */
{
  short i;
  char *ptr, *ptre, *ptr2;
  for (i = 1; i <= dim; i++) {
    ptr = a[i];
    ptre = ptr + dim;
    ptr2 = b[i];
    while (++ptr <= ptre)
      *(++ptr2) = *ptr;
  }
}

int ncopy(int n, char ** a, char ** b)
/* Copies n times matrix a to matrix b. Externals: prime,dim. */
{
  char *ptr, *ptre, *ptr2;
  short x, i;
  for (i = 1; i <= dim; i++) {
    ptr = a[i];
    ptre = ptr + dim;
    ptr2 = b[i];
    while (++ptr <= ptre) {
      x = n * *ptr;
      *(++ptr2) = x % prime;
    }
  }
}

int sum(int n, char ** a, char ** b)
/* Adds n times matrix a to matrix b. Externals: prime,dim. */
{
  char *ptr, *ptre, *ptr2;
  int   x;
  short i;
  for (i = 1; i <= dim; i++) {
    ptr = a[i];
    ptre = ptr + dim;
    ptr2 = b[i];
    while (++ptr <= ptre) {
      x = *(++ptr2) + n * *ptr;
      *ptr2 = x % prime;
    }
  }
}

int im(char * v, char * w, char ** a)
/* The image w of vector v under matrix a is computed.
   Externals: dim,prime.
*/
{
  char *p, *q;
  short i, j;
  int   x;
  p = w;
  q = w + dim;
  j = 0;
  while (++p <= q) {
    j++;
    x = 0;
    for (i = 1; i <= dim; i++)
      x += (v[i] * a[i][j]);
    *p = x % prime;
  }
}

int prod(char ** a, char ** b, char ** c)
/* The product of matrices a and b is computed and stored in c.
   Warning:  a and c may be equal, but not b and c.
   Externals: prime,dim,cvec.
*/
{
  char *ptr, *ptre, *ptr2;
  short i;
  for (i = 1; i <= dim; i++) {
    im(a[i], cvec, b);
    ptr = cvec;
    ptre = ptr + dim;
    ptr2 = a[i];
    while (++ptr <= ptre)
      *(++ptr2) = *ptr;
  }
}

int readmat(char ** a)
/* Matrix a is read from input ip.
   Externals: dim,ip.
*/
{
  char *ar, *are;
  short rno, i;
  for (i = 1; i <= dim; i++) {
    ar = a[i];
    are = ar + dim;
    while (++ar <= are) {
      fscanf(ip, "%hd", &rno);
      *ar = rno;
    }
  }
}

int rvecsum(int n, char * v)
/* n times vector read from ip is added to v. Externals: prime,dim,ip. */
{
  char *ptr, *ptre;
  int   x;
  short rno;
  ptr = v;
  ptre = ptr + dim;
  while (++ptr <= ptre) {
    fscanf(ip, "%hd", &rno);
    x = *ptr + n * rno;
    *ptr = x % prime;
  }
}

int printmat(char ** a)
/* Matrix a is output to op.
   Externals: dim,op.
*/
{
  char *ar, *are;
  short i;
  for (i = 1; i <= dim; i++) {
    ar = a[i];
    are = ar + dim;
    if (prime < 10)
      while (++ar <= are)
        fprintf(op, "%2d", *ar);
    else if (prime < 100)
      while (++ar <= are)
        fprintf(op, "%3d", *ar);
    else
      while (++ar <= are)
        fprintf(op, "%4d", *ar);
    fprintf(op, "\n");
  }
  fprintf(op, "\n");
}

int null(char ** a, char ** b)
/* Generators of the null space of matrix a are output to op.
   Matrices a and b are destroyed by this process.
   Rank of nullspace is returned.
   If aut is set, then we exit as soon as nullity exceeds maxnull.
   Externals:  prime,dim,op,aut,maxnull.
*/
{
  char *ai, *br, *ptr, *ptre, *ptr2, *ptr3, *ptr4, x, y;
  int   i, j, k, z, nlty;
  for (i = 1; i <= dim; i++) {
    br = b[i];
    ptr = br;
    ptre = ptr + dim;
    while (++ptr <= ptre)
      *ptr = 0;
    br[i] = 1;
  }
  nlty = 0;
  for (i = 1; i <= dim; i++) {
    ai = a[i];
    for (j = 1; j <= dim; j++)
      if ((x = ai[j]) != 0)
        break;
    if (j <= dim) {
      y = pinv[x];
      ptr = ai;
      ptre = ptr + dim;
      ptr2 = b[i];
      while (++ptr <= ptre) {
        z = *ptr * y;
        *ptr = z % prime;
        z = *(++ptr2) * y;
        *ptr2 = z % prime;
      }
      for (k = i + 1; k <= dim; k++) {
        ptr2 = a[k];
        if ((x = ptr2[j]) != 0) {
          x = prime - x;
          ptr = ai;
          ptre = ptr + dim;
          ptr3 = b[k];
          ptr4 = b[i];
          while (++ptr <= ptre) {
            z = *(++ptr2) + x * *ptr;
            *ptr2 = z % prime;
            z = *(++ptr3) + x * *(++ptr4);
            *ptr3 = z % prime;
          }
        }
      }
    }
    else {
      nlty++;
      if (aut && nlty > maxnull) {
        printf("Nullity>%d\n", maxnull);
        fflush(stdout);
        return (nlty);
      }
      ptr = b[i];
      ptre = ptr + dim;
      if (prime < 10)
        while (++ptr <= ptre)
          fprintf(op, "%2d", *ptr);
      else if (prime < 100)
        while (++ptr <= ptre)
          fprintf(op, "%3d", *ptr);
      else
        while (++ptr <= ptre)
          fprintf(op, "%4d", *ptr);
      fprintf(op, "\n");
    }
  }
  printf("Null space has dimension %d.\n", nlty);
  fflush(stdout);
  return (nlty);
}

int null1(char ** a, char ** b)
/* One generator of the null space of matrix a is output to op.
   Matrices a and b are destroyed by this process.
   Externals:  prime,dim,op.
*/
{
  char *ai, *br, *ptr, *ptre, *ptr2, *ptr3, *ptr4, x, y;
  int   i, j, k, z;
  for (i = 1; i <= dim; i++) {
    br = b[i];
    ptr = br;
    ptre = ptr + dim;
    while (++ptr <= ptre)
      *ptr = 0;
    br[i] = 1;
  }
  for (i = 1; i <= dim; i++) {
    ai = a[i];
    for (j = 1; j <= dim; j++)
      if ((x = ai[j]) != 0)
        break;
    if (j <= dim) {
      y = pinv[x];
      ptr = ai;
      ptre = ptr + dim;
      ptr2 = b[i];
      while (++ptr <= ptre) {
        z = *ptr * y;
        *ptr = z % prime;
        z = *(++ptr2) * y;
        *ptr2 = z % prime;
      }
      for (k = i + 1; k <= dim; k++) {
        ptr2 = a[k];
        if ((x = ptr2[j]) != 0) {
          x = prime - x;
          ptr = ai;
          ptre = ptr + dim;
          ptr3 = b[k];
          ptr4 = b[i];
          while (++ptr <= ptre) {
            z = *(++ptr2) + x * *ptr;
            *ptr2 = z % prime;
            z = *(++ptr3) + x * *(++ptr4);
            *ptr3 = z % prime;
          }
        }
      }
    }
    else {
      ptr = b[i];
      ptre = ptr + dim;
      if (prime < 10)
        while (++ptr <= ptre)
          fprintf(op, "%2d", *ptr);
      else if (prime < 100)
        while (++ptr <= ptre)
          fprintf(op, "%3d", *ptr);
      else
        while (++ptr <= ptre)
          fprintf(op, "%4d", *ptr);
      fprintf(op, "\n");
      return (1);
    }
  }
  printf("Null space is trivial.\n");
  return (0);
}

int spgen(char ** a, int n)
/* The space generated by the images of vector a[1] under mats 1-n is
   computed. Normed gens of this space are put in initial rows of matrix a,
   and then the remaining rows are added to make a nonsingular. Externals:
   prime,dim,svec,mat.
*/
{
  int   i, j, k, z, rno, spdim;
  char *ar, *cf, *ptr, *ptre, *ptr2, l, x;
  ptr = a[1];
  ptre = ptr + dim;
  while (++ptr <= ptre)
    if (*ptr != 0)
      break;
  if (ptr > ptre) {
    printf("Vector is zero.\n");
    return (0);
  }
  svec[1] = ptr - a[1];
  x = pinv[*ptr];
  ptr--;
  while (++ptr <= ptre) {
    z = x * *ptr;
    *ptr = z % prime;
  }
  spdim = 1;
  for (i = 1; i <= spdim; i++)
    for (j = 1; j <= n; j++) {
      ar = a[spdim + 1];
      im(a[i], ar, mat[j]);
      for (k = 1; k <= spdim; k++) {
        if ((l = ar[svec[k]]) != 0) {
          l = prime - l;
          ptr = ar;
          ptre = ptr + dim;
          ptr2 = a[k];
          while (++ptr <= ptre) {
            z = *ptr + l * *(++ptr2);
            *ptr = z % prime;
          }
        }
      }
      ptr = ar;
      ptre = ptr + dim;
      while (++ptr <= ptre)
        if (*ptr != 0)
          break;
      if (ptr <= ptre) {
        spdim++;
        svec[spdim] = ptr - ar;
        x = pinv[*ptr];
        ptr--;
        while (++ptr <= ptre) {
          z = x * *ptr;
          *ptr = z % prime;
        }
        if (spdim == dim) {
          printf("Whole space is generated.\n");
          fflush(stdout);
          return (dim);
        }
      }
    }
  printf("Subspace generated has dimension %d.\n", spdim);
  fflush(stdout);
  cf = a[dim];
  rno = spdim;
  ptr = cf;
  ptre = ptr + dim;
  while (++ptr <= ptre)
    *ptr = 0;
  for (i = 1; i <= spdim; i++)
    cf[svec[i]] = 1;
  for (i = 1; i <= dim; i++)
    if (cf[i] == 0) {
      rno++;
      svec[rno] = i;
      ar = a[rno];
      ptr = ar;
      ptre = ar + dim;
      while (++ptr <= ptre)
        *ptr = 0;
      ar[i] = 1;
      if (rno == dim)
        break;
    }
  return (spdim);
}

int opnmat(char ** a, int n, int tdim, int fop)
/* Output of matrices 1-n using basis in mat a.
   tdim=dim of output matrices. Output coeffs begin at fop.
   Externals: prime,dim,mat,cvec,svec,op.
*/
{
  int  i, j, k, z;
  char c, *ptr, *ptre, *ptr2;
  for (i = 1; i <= n; i++)
    for (j = fop; j <= tdim; j++) {
      ptre = cvec + dim;
      im(a[j], cvec, mat[i]);
      for (k = 1; k <= tdim; k++) {
        c = cvec[svec[k]];
        if (k >= fop) {
          if (prime < 10)
            fprintf(op, "%2d", c);
          else if (prime < 100)
            fprintf(op, "%3d", c);
          else
            fprintf(op, "%4d", c);
        }
        if (c > 0) {
          c = prime - c;
          ptr = cvec;
          ptr2 = a[k];
          while (++ptr <= ptre) {
            z = *ptr + c * *(++ptr2);
            *ptr = z % prime;
          }
        }
      }
      fprintf(op, "\n");
    }
}
