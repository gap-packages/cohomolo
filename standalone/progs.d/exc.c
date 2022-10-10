#include "defs.h"

extern char  mult, outft[];
extern short rwd[], ***scoeff[], *cst, *cend, *ocst, ***coeff, wd1[], wd2[],
    wd3[], wd4[];
extern int     rsp;
extern short **simcos[], pinv[], gno[], base[], lorb[], *svptr[], cord[],
    cp[], bno, invg[], **mat[], mwl2, mwdl, mlwdl, maxcos, dim, prime, *spv,
    **spm, ng, rwl, nb;
extern FILE *ip, *op;

int gc(void)
/* Garbage collection of coefficients by writing and rereading */
{
  short i, j;
  short l, *p, *q;
  op = fopen(outft, "w");
  for (i = 0; i < ng; i++)
    for (j = 1; j <= maxcos; j++)
      if ((p = coeff[i][j]) != 0) {
        l = *p;
        q = p + l;
        while (p <= q) {
          fprintf(op, " %d", *p);
          p++;
        }
        l = *p;
        q = p + l;
        while (p <= q) {
          fprintf(op, " %d", *p);
          p++;
        }
      }
  fclose(op);
  ip = fopen(outft, "r");
  cst = ocst;
  for (i = 0; i < ng; i++)
    for (j = 1; j <= maxcos; j++)
      if (coeff[i][j] != 0) {
        coeff[i][j] = cst;
        fscanf(ip, "%hd", cst);
        p = cst;
        q = p + *p;
        while (++p <= q)
          fscanf(ip, "%hd", p);
        cst = q + 1;
        fscanf(ip, "%hd", cst);
        p = cst;
        q = p + *p;
        while (++p <= q)
          fscanf(ip, "%hd", p);
        cst = q + 1;
      }
  fclose(ip);
  unlink("gctemp");
  rsp = cend + 1 - cst;
  printf("Performed garbage collection. rsp=%d.\n", rsp);
  return (0);
}

int comp(short * a, short ** dp)
/* Compress string + vector in "a" into string in *dp
   Vectors always occupy the last dim entries of words.
*/
{
  short *a1, ol, nl, nl2, *p, *q, *r, *ae;
  a1 = a + mwdl;
  ae = a1 + dim;
  p = a1;
  nl2 = 0;
  while (++p <= ae)
    if (*p != 0)
      nl2 += 2;
  ol = (*dp == 0) ? 0 : **dp + *(*dp + **dp + 1);
  nl = nl2 + *a;
  if (nl == 0) {
    *dp = 0;
    return (0);
  }
  if (nl > ol) {
    nl += 2;
    if (rsp < nl)
      gc();
    if (rsp < nl) {
      fprintf(stderr, "Out of cspace. Increase CSPACE.\n");
      return (-1);
    }
    *dp = cst;
    cst += nl;
    rsp -= nl;
  }
  p = *dp;
  q = a;
  r = q + *a;
  while (q <= r) {
    *p = *q;
    p++;
    q++;
  }
  *p = nl2;
  q = a1;
  while (++q <= ae)
    if (*q != 0) {
      *(++p) = q - a1;
      *(++p) = *q;
    }
  return (0);
}

int compb(short * a, short ** dp)
/* Similar, but rewrites word at end of short space */
{
  short *a1, nl, nl2, *p, *q, *r, *ae;
  a1 = a + mwdl;
  ae = a1 + dim;
  p = a1;
  nl2 = 0;
  while (++p <= ae)
    if (*p != 0)
      nl2 += 2;
  nl = nl2 + *a;
  if (nl == 0) {
    *dp = 0;
    return (0);
  }
  nl += 2;
  if (rsp < nl)
    gc();
  if (rsp < nl) {
    fprintf(stderr, "Out of cspace(b). Increase CSPACE.\n");
    return (-1);
  }
  cend -= nl;
  rsp -= nl;
  *dp = cend + 1;
  p = *dp;
  q = a;
  r = q + *a;
  while (q <= r) {
    *p = *q;
    p++;
    q++;
  }
  *p = nl2;
  q = a1;
  while (++q <= ae)
    if (*q != 0) {
      *(++p) = q - a1;
      *(++p) = *q;
    }
  return (0);
}

int expand(short * a, short * b)
/* expand string b to string + vector a */
{
  short *a1, *ae, *p, *q, *r;
  a1 = a + mwdl;
  ae = a1 + dim;
  p = a1;
  while (++p <= ae)
    *p = 0;
  if (b == 0) {
    *a = 0;
    return (0);
  }
  *a = *b;
  p = a;
  q = a + *a;
  r = b;
  while (++p <= q)
    *p = *(++r);
  r++;
  q = r + *r;
  while (++r < q) {
    a1[*r] = *(r + 1);
    r++;
  }
  return (0);
}

int action(short * a, short * b)
/* Replace vector a by its image under action by matrices of gens in word b
 */
{
  short  z, *be, c, *p, *ae;
  short *spve, *s, *v;
  z = 1;
  p = a;
  ae = a + dim;
  while (++p <= ae)
    if (*p != 0) {
      z = 0;
      break;
    }
  if (z)
    return (0);
  spve = spv + dim;
  be = b + *b;
  while (++b <= be) {
    s = spv;
    while (++s <= spve)
      *s = 0;
    p = a;
    while (++p <= ae)
      if ((c = *p) != 0) {
        v = mat[*b][p - a];
        s = spv;
        while (++s <= spve) {
          *s += (c * *(++v));
          *s %= prime;
        }
      }
    p = a;
    s = spv;
    while (++p <= ae)
      *p = *(++s);
  }
  return (0);
}

int concat(short * a, short * b)
/* Concatenate strings+vectors "a" and "b" */
{
  short *a1, *bb, *be, c, *p, *q, *r, canct, min;
  if (b == 0)
    return (0);
  a1 = a + mwdl;
  bb = b + *b + 1;
  be = bb + *bb;
  if (mult) {
    while (++bb < be) {
      c = *bb;
      bb++;
      a1[c] += *bb;
      a1[c] %= cord[c];
    }
  }
  else {
    action(a1, b);
    while (++bb < be) {
      c = *bb;
      bb++;
      a1[c] += *bb;
      a1[c] %= prime;
    }
  }
  min = (*a <= *b) ? *a : *b;
  canct = 0;
  p = a + *a;
  q = b + 1;
  while (canct < min && *p == invg[*q]) {
    canct++;
    p--;
    q++;
  }
  r = b + *b;
  while (q <= r)
    *(++p) = *(q++);
  *a += (*b - 2 * canct);
  if (*a > mwl2)
    if (reduce(a) == -1)
      return (-1);
  return (0);
}

int concatl(short * b)
/* Similar, but concatenate b to the long word rwd */
{
  short *a1, *bb, *be, c, *p, *q, *r, canct, min;
  if (b == 0)
    return (0);
  a1 = rwd + mlwdl;
  bb = b + *b + 1;
  be = bb + *bb;
  if (mult) {
    while (++bb < be) {
      c = *bb;
      bb++;
      a1[c] += *bb;
      a1[c] %= cord[c];
    }
  }
  else {
    action(a1, b);
    while (++bb < be) {
      c = *bb;
      bb++;
      a1[c] += *bb;
      a1[c] %= prime;
    }
  }
  min = (rwl <= *b) ? rwl : *b;
  canct = 0;
  p = rwd + rwl;
  q = b + 1;
  while (canct < min && *p == invg[*q]) {
    canct++;
    p--;
    q++;
  }
  r = b + *b;
  rwl += (*b - 2 * canct);
  if (rwl > mlwdl) {
    fprintf(stderr, "Word too long. Increase MLWDL.\n");
    return (-1);
  }
  while (q <= r)
    *(++p) = *(q++);
  return (0);
}

int invwd(short * a, short * b)
/* Invert the word+vector a into b */
{
  short *a1, *b1, *be, *p, *q;
  a1 = a + mwdl;
  p = a1;
  b1 = b + mwdl;
  be = b1 + dim;
  q = b1;
  while (++q <= be)
    *q = 0;
  q = b1;
  if (mult) {
    while (++q <= be)
      if (*(++p) != 0)
        *q = cord[p - a1] - *p;
  }
  else
    while (++q <= be)
      if (*(++p) != 0)
        *q = prime - *p;
  *b = *a;
  p = a + *a;
  q = b + 1;
  while (p > a) {
    *q = invg[*p];
    p--;
    q++;
  }
  if (mult == 0)
    action(b1, b);
  return (0);
}

int ainvb(short * a, short * b, short * c)
/* Concatenate word+vector "a" with inverse of word+vector "b" and write into
 * c
 */
{
  short min, canct, *p, *q, *r, *c1, *ce, *a1;
  min = (*a <= *b) ? *a : *b;
  p = a + 1;
  q = b + 1;
  canct = 0;
  while (canct < min && *p == *q) {
    canct++;
    p++;
    q++;
  }
  p = a + *a + 1;
  q = a + canct + 1;
  r = c;
  *r = *a + *b - 2 * canct;
  while (--p >= q)
    *(++r) = invg[*p];
  p = b + canct;
  q = b + *b;
  while (++p <= q)
    *(++r) = *p;
  a1 = a + mwdl;
  c1 = c + mwdl;
  ce = c1 + dim;
  p = c1;
  while (++p <= ce)
    *p = 0;
  p = c1;
  q = a1;
  r = b + mwdl;
  if (mult) {
    while (++p <= ce)
      if (*(++q) != 0)
        *p = cord[p - c1] - *q;
    p = c1;
    while (++p <= ce)
      if (*(++r) != 0) {
        *p += *r;
        *p %= cord[p - c1];
      }
  }
  else {
    while (++p <= ce)
      if (*(++q) != 0)
        *p = prime - *q;
    action(c1, c);
    p = c1;
    while (++p <= ce)
      if (*(++r) != 0) {
        *p += *r;
        *p %= prime;
      }
  }
  if (*c > mwl2)
    if (reduce(c) == -1)
      return (-1);
  return (0);
}

int reduce(short * a)
/* Attempt to reduce length of word+vector "a", using coefficients computed
   from higher values of bno
*/
{
  short *q, *r, *a1, *rwd1, *rwde;
  short  ol, nl, cos, i, j, *p, *s, *cp1, *cpe;
  ol = *a;
  *cp = ol;
  p = cp;
  q = a;
  r = a + ol;
  cp1 = cp + mlwdl;
  cpe = cp1 + dim;
  while (++q <= r)
    *(++p) = *q;
  p = cp1;
  while (++p <= cpe)
    *p = 0;
  for (i = bno + 1; i <= nb; i++)
    addsv(image(base[i]), svptr[i]);
  nl = *cp - ol;
  if (nl < ol) {
    *a = nl;
    p = cp + *cp + 1;
    q = a;
    r = a + nl;
    while (++q <= r)
      *q = invg[*(--p)];
    rwd1 = rwd + mlwdl;
    rwde = rwd1 + dim;
    for (i = bno + 1; i <= nb; i++)
      if (lorb[i] > 1) {
        q = rwd1;
        while (++q <= rwde)
          *q = 0;
        cos = 1;
        rwl = 0;
        p = cp;
        s = p + *cp;
        while (++p <= s) {
          if (concatl(scoeff[i][*p][cos]) == -1)
            return (-1);
          cos = simcos[i][*p][cos];
        }
        if (cos != 1) {
          fprintf(stderr, "Reduction error.\n");
          return (-1);
        }
        *cp = rwl;
        p = cp;
        s = p + rwl;
        q = rwd;
        while (++p <= s)
          *p = *(++q);
        q = rwd1;
        r = rwde;
        p = cp1;
        if (mult)
          while (++q <= r) {
            *(++p) += *q;
            j = cord[q - rwd1];
            *p %= j;
          }
        else
          while (++q <= r) {
            *(++p) += *q;
            *p %= prime;
          }
      }
    p = cp1;
    r = rwd1;
    while (++p <= cpe)
      *(++r) = *p;
    a1 = a + mwdl;
    r = rwd1;
    q = a1;
    if (mult) {
      while (++r <= rwde) {
        *(++q) += *r;
        *q %= cord[q - a1];
      }
    }
    else {
      action(rwd1, a);
      while (++r <= rwde) {
        *(++q) += *r;
        *q %= prime;
      }
    }
  }
  if (*a > mwl2) {
    fprintf(stderr, "Reduction unsuccessful.\n");
    return (-1);
  }
  return (0);
}

int concheck(void)
/* Consistency check */
{
  short *p, *q, ok;
  if (ainvb(wd1, wd2, wd3) == -1)
    return (-1);
  if (reduce(wd3) == -1)
    return (-1);
  ok = *wd3 == 0;
  if (ok) {
    p = wd3 + mwdl;
    q = p + dim;
    while (++p <= q)
      if (*p != 0) {
        ok = 0;
        break;
      }
  }
  if (ok == 0) {
    fprintf(stderr, "Inconsistent relation.\n");
    return (-1);
  }
  return (0);
}

int setpinv(void)
{
  short i, j;
  for (i = 0; i < prime; i++)
    pinv[i] = 0;
  for (i = 1; i < prime; i++)
    if (pinv[i] == 0)
      for (j = 1; j < prime; j++)
        if (i * j % prime == 1) {
          pinv[i] = j;
          pinv[j] = i;
          break;
        }
  return (0);
}
