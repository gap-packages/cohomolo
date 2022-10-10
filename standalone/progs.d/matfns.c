#include "defs.h"

extern short prime, dim, *spv, **spm, **mat[], pinv[];
extern FILE *ip, *op;

void trans(short ** a, short ** b)
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

void im(short * v, short * w, short ** a)
/* The image w of vector v under matrix a is computed.
   Externals: dim,prime.
*/
{
  short *p, *q, i, j;
  int    sum;
  p = w;
  q = w + dim;
  j = 0;
  while (++p <= q) {
    j++;
    sum = 0;
    for (i = 1; i <= dim; i++)
      sum += (v[i] * a[i][j]);
    *p = sum % prime;
  }
}

int comm(short * v, short * w, short ** a)
/* The image of v under a minus v is computed into w. The position of the
   first nonzero entry in w is returned, or 0 if w=0.
   Externals:  dim,prime.
*/
{
  short *p, *q, i, j, k;
  int    sum;
  p = w;
  q = w + dim;
  j = 0;
  k = 0;
  while (++p <= q) {
    j++;
    sum = -v[j];
    for (i = 1; i <= dim; i++)
      sum += (v[i] * a[i][j]);
    *p = sum % prime;
    if (*p < 0)
      *p += prime;
    if (k == 0 && *p != 0)
      k = p - w;
  }
  return (k);
}

void prod(short * cm, short ** a)
/* The product of matrices mat[cm[1]],...mat[cm[*cm]] is computed into a.
   Externals:  mat,dim.
*/
{
  short i, j, l, m, *v, *w, *x;
  ;
  l = *cm;
  m = l % 2;
  if (l == 0)
    for (i = 1; i <= dim; i++) {
      v = a[i];
      w = v + dim;
      while (++v <= w)
        *v = 0;
      a[i][i] = 1;
    }
  else if (l == 1)
    for (i = 1; i <= dim; i++) {
      v = a[i];
      w = v + dim;
      x = mat[cm[1]][i];
      while (++v <= w) {
        x++;
        *v = *x;
      }
    }
  else
    for (i = 1; i <= dim; i++) {
      v = a[i];
      if (m == 0)
        im(mat[cm[1]][i], v, mat[cm[2]]);
      else {
        im(mat[cm[1]][i], spv, mat[cm[2]]);
        im(spv, v, mat[cm[3]]);
      }
      for (j = 3 + m; j <= l; j += 2) {
        im(v, spv, mat[cm[j]]);
        im(spv, v, mat[cm[j + 1]]);
      }
    }
}

int inv(short ** a, short ** b)
/* Matrix a is inverted into b. -1 is returned if a is singular, otherwise 0.
   Externals: spm,prime,dim.
*/
{
  short i, j, fac, *ar, *are, *br, *spmr, *spme, *spmv;
  int   sum;
  for (i = 1; i <= dim; i++) {
    ar = a[i];
    br = b[i];
    spmr = spm[i];
    are = ar + dim;
    while (++ar <= are) {
      spmr++;
      br++;
      *spmr = *ar;
      *br = 0;
    }
    b[i][i] = 1;
  }
  for (i = 1; i <= dim; i++) {
    spmr = spm[i];
    spme = spmr + dim;
    br = b[i];
    if (pinv[spmr[i]] == 0) {
      for (j = i + 1; j <= dim; j++)
        if (pinv[spm[j][i]] != 0)
          break;
      if (j > dim) {
        fprintf(stderr, "Matrix is singular.\n");
        return (-1);
      }
      spmr = spm[j];
      spm[j] = spm[i];
      spm[i] = spmr;
      spme = spmr + dim;
      br = b[j];
      b[j] = b[i];
      b[i] = br;
    }
    fac = pinv[spmr[i]];
    spmr++;
    br++;
    while (spmr <= spme) {
      sum = *spmr * fac;
      *spmr = sum % prime;
      sum = *br * fac;
      *br = sum % prime;
      spmr++;
      br++;
    }
    for (j = 1; j <= dim; j++) {
      if ((fac = spm[j][i]) == 0 || j == i)
        continue;
      spmr = spm[i] + 1;
      spmv = spm[j] + 1;
      br = b[i] + 1;
      ar = b[j] + 1;
      while (spmr <= spme) {
        sum = *spmv - (*spmr * fac);
        *spmv = sum % prime;
        if (*spmv < 0)
          *spmv += prime;
        sum = *ar - (*br * fac);
        *ar = sum % prime;
        if (*ar < 0)
          *ar += prime;
        spmr++;
        spmv++;
        br++;
        ar++;
      }
    }
  }
  return (0);
}

void readmat(short ** a)
/* Matrix a is read from input ip.
   Externals: dim,ip,prime.
*/
{
  short i, *ar, *are;
  for (i = 1; i <= dim; i++) {
    ar = a[i];
    are = ar + dim;
    while (++ar <= are) {
      fscanf(ip, "%hd", ar);
      while (*ar < 0)
        *ar += prime;
      *ar %= prime;
    }
  }
}

void printmat(short ** a)
/* Matrix a is output to op.
   Externals: dim,op.
*/
{
  short i, *ar, *are;
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

int cbdef(
    int gb, int ge, int cbno, short * d1, short * d2, short * wt, short * acl)
/* It is assumed that matrices mat[gb],...mat[ge] generate a p-group,
   and that the action of this group on the vector space has already been
   triangularized. (p=prime).
    A basis change matrix is computed, as mat[cbno] (rows=new basis
   elts in terms of old), to make this action on the space suitable for a PCP.
   The weights of the basis elts 1,..,dim of the space under the action are
   returned as wt[1],...wt[dim], and the class as *acl. The definitions of the
   basis elts of weight>1 will be given by the last term of [d1[i],d2[i]],
  where d1[i] is a basis elt, and d2[i] a matrix no. 1 is returned if
  mat[cbno] is the identity, otherwise 0. Externals: mat,spv,prime,dim.
*/
{
  short i, j, k, l, fac, *ar, *are, *p, **cbm;
  int   sum;
  char  id;
  cbm = mat[cbno];
  for (i = 1; i <= dim; i++) {
    ar = cbm[i];
    are = ar + dim;
    p = ar;
    while (++p <= are)
      *p = 0;
    ar[i] = 1;
    d1[i] = 0;
    d2[i] = 0;
    wt[i] = 1;
  }
  id = 1;
  *acl = 1;
restart:
  for (i = 1; i <= dim; i++)
    for (j = gb; j <= ge; j++) {
      k = comm(cbm[i], spv, mat[j]);
      while (k > 0)
        if (wt[k] <= wt[i] || (d1[k] == 0 && wt[k] == wt[i] + 1)) {
          wt[k] = wt[i] + 1;
          if (wt[k] > *acl)
            *acl = wt[k];
          if (id && spv[k] != 1)
            id = 0;
          if (id)
            for (l = k + 1; l <= dim; l++)
              if (spv[l] > 0)
                id = 0;
          if (id == 0) {
            ar = cbm[k] + k - 1;
            are = cbm[k] + dim;
            p = spv + k - 1;
            while (++ar <= are)
              *ar = *(++p);
            if (d1[k] > 0) {
              for (l = 1; l <= dim; l++) {
                d1[l] = 0;
                d2[l] = 0;
              }
              goto restart;
            }
          }
          d1[k] = i;
          d2[k] = j + 1 - gb;
          k = 0;
        }
        else {
          ar = cbm[k] + k;
          are = ar + dim - k;
          sum = spv[k] * pinv[*ar];
          fac = sum % prime;
          fac = prime - fac;
          p = spv + k;
          k = 0;
          while (ar <= are) {
            sum = *p + (fac * *ar);
            *p = sum % prime;
            if (k == 0 && *p != 0)
              k = p - spv;
            p++;
            ar++;
          }
        }
    }
  return (id);
}
