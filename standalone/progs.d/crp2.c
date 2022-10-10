#include "defs.h"

extern char   mult, inf0[], inf1[], inf2[], outf[], outft[], inf3[];
extern short *cst, **cpst, ***cdpst, csp[], *cpsp[], **cdpsp[], ***coeff[];
extern int    space, cptrsp, cspace;
extern short  sp[], **mat[], *psp[], **imcos[], **cpco[], lorb[], pinv[], wdl,
    ptrsp, cdptrsp, mpr, marg, *spst, **pspst, npt, nb, nph, nph2, rno, orno,
    coh_index, *invg;
short *wd1, *wd2;
short *gno, *cord, cl, dim, prime, **vec, **spm, *spv, wv, *val;
FILE * ipr, *op, *ipx, *ip;

int crprog2(void)
/* In the second half, the values of the words in H in the module or
   multiplier are computed, using the coefficients computed in extprun.
*/
{
  short  i, j, k, l, m, n, rct, inct, x, l1, l2, *ps, *qs;
  short *p, *q, *r, *v2, *ex, c;
  int    prod;
  printf("Beginning crprog2.\n");
  spst = sp + nph2;
  pspst = psp;
  strcpy(inf0, inf2);
  if (mult == 0) {
    strcat(inf2, "mat");
    if ((ip = fopen(inf2, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf2);
      return (-1);
    }
    fscanf(ip, "%hd%hd%hd", &prime, &dim, &k);
    if (prime > mpr) {
      fprintf(stderr, "Prime too big. Increase MPR.\n");
      return (-1);
    }
    if (k != nph) {
      fprintf(stderr, "No of mats wrong.\n");
      return (-1);
    }
    setpinv();
    vec = psp;
    n = (nph2 + 1) * dim;
    pspst += n;
    if (pspst - psp >= ptrsp) {
      fprintf(stderr, "Out of ptr space. Increase PTRSP.\n");
      return (-1);
    }
    for (i = 0; i < n; i++)
      vec[i] = spst - 1 + i * dim;
    spst += (dim * n);
    if (spst - sp > space) {
      fprintf(stderr, "Out of space. Increase SPACE.\n");
      return (-1);
    }
    for (i = 0; i <= nph2; i++)
      mat[i] = vec - 1 + i * dim;
    spm = mat[nph2];
    spv = spm[1];
    for (i = 0; i < nph2; i += 2) {
      readmat(mat[i]);
      inv(mat[i], mat[i + 1]);
    }
    fclose(ip);
  }

  /* Now we read in the output of extprun */
  wd1 = csp - 1;
  wd2 = wd1 + wdl;
  cst = csp + 2 * wdl;
  cpst = cpsp;
  cdpst = cdpsp;
  strcpy(inf2, inf0);
  strcat(inf2, ".ep");
  if ((ip = fopen(inf2, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf2);
    return (-1);
  }
  fscanf(ip, "%hd%hd", &i, &j);
  if (i != nb || (mult == 0 && j != dim)) {
    fprintf(stderr, "%s has line 1 wrong.\n", inf2);
    return (-1);
  }
  if (mult) {
    dim = j;
    cord = spst - 1;
    spst += dim;
  }
  gno = spst - 1;
  spst += nb;
  val = spst - 1;
  spst += dim;
  for (i = 1; i <= nb; i++) {
    fscanf(ip, "%hd", &j);
    if (j != lorb[i]) {
      fprintf(stderr, "lorb wrong in %s.\n", inf2);
      return (-1);
    }
  }
  if (mult)
    for (i = 1; i <= dim; i++)
      fscanf(ip, "%hd", cord + i);
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%hd", gno + i);
  for (i = 1; i <= nb; i++) {
    n = gno[i];
    if ((l = lorb[i]) > 1) {
      if (pspst + n - psp > ptrsp) {
        fprintf(stderr, "Out of ptrsp. Increase PTRSP.\n");
        return (-1);
      }
      if (cdpst + n - cdpsp > cdptrsp) {
        fprintf(stderr, "Out of cdsp. Increase CDPTRSP.\n");
        return (-1);
      }
      prod = n * l;
      if (spst + prod - sp > space) {
        fprintf(stderr, "Out of space. Increase SPACE.\n");
        return (-1);
      }
      if (cpst + prod - cpsp > cptrsp) {
        fprintf(stderr, "Out of cptrsp. Increase CPTRSP.\n");
        return (-1);
      }
      imcos[i] = pspst;
      pspst += n;
      coeff[i] = cdpst;
      cdpst += n;
      for (j = 0; j < n; j++) {
        imcos[i][j] = spst - 1;
        for (k = 1; k <= l; k++) {
          fscanf(ip, "%hd", spst);
          spst++;
        }
        coeff[i][j] = cpst - 1;
        for (k = 1; k <= l; k++) {
          fscanf(ip, "%hd", &l1);
          if (l1 == 0) {
            fscanf(ip, "%hd", &l2);
            if (l2 == 0)
              *(cpst++) = 0;
            else {
              *(cpst++) = cst;
              *(cst++) = 0;
              *(cst++) = l2;
              for (m = 1; m <= l2; m++) {
                fscanf(ip, "%hd", &x);
                *(cst++) = x;
              }
            }
          }
          else {
            *(cpst++) = cst;
            *(cst++) = l1;
            for (m = 1; m <= l1; m++) {
              fscanf(ip, "%hd", &x);
              *(cst++) = x;
            }
            fscanf(ip, "%hd", &l2);
            *(cst++) = l2;
            for (m = 1; m <= l2; m++) {
              fscanf(ip, "%hd", &x);
              *(cst++) = x;
            }
          }
        }
        if (cst + marg - csp > cspace) {
          fprintf(stderr, "Running out of cspace. Increase CSP.\n");
          return (-1);
        }
      }
    }
  }
  fclose(ip);
  /* End of input. Now computation can begin */

  op = fopen(outf, "w");
  ipr = fopen(inf1, "r");
  fscanf(ipr, "%hd", &i);
  fprintf(op, "%4d%4d\n", nb, dim);
  for (i = 1; i <= nb; i++) {
    fscanf(ipr, "%hd", &j);
    fprintf(op, " %4d", j);
  }
  while (getc(ipr) != '\n')
    ;
  fprintf(op, "\n");
  if (mult) {
    for (i = 1; i <= dim; i++)
      fprintf(op, "%4d", cord[i]);
    fprintf(op, "\n");
  }
  wv = wdl - dim - 1;
  ipx = fopen(outft, "r");

  for (rct = 1; rct <= rno; rct++) {
    printf("Scanning relation no %d\n", rct);
    if (rct == orno + 1) {
      while ((c = getc(ipr)) != '\n')
        putc(c, op);
      putc('\n', op);
    }
    while ((c = getc(ipr)) != '\n')
      putc(c, op);
    putc('\n', op);
    if (mult == 0) {
      if ((ip = fopen(inf3, "r")) == 0) {
        fprintf(stderr, "Cannot open %s.\n", inf3);
        return (-1);
      }
      fscanf(ip, "%hd%hd%hd", &i, &j, &k);
      if (i != prime || j != dim || k != (coh_index - 1)) {
        fprintf(stderr, "%s has line 1 wrong.\n", inf3);
        return (-1);
      }
    }
    for (j = 1; j <= dim; j++)
      val[j] = 0;
    for (inct = 1; inct <= coh_index; inct++) {
      fscanf(ipx, "%hd", &l2);
      for (k = 1; k <= l2; k++) {
        fscanf(ipx, "%hd", &l);
        wd2[k] = l;
      }
      p = wd1 + wv;
      q = wd2 + wv;
      r = p + dim;
      while (++p <= r) {
        *p = 0;
        *(++q) = 0;
      }
      for (k = 1; k <= nb; k++)
        if (lorb[k] > 1) {
          if (scan(k, &l2) == -1)
            return (-1);
          ex = wd1;
          wd1 = wd2;
          wd2 = ex;
        }
      v2 = wd2 + wv;
      ps = val;
      if (inct > 1 && mult == 0) {
        readmat(spm);
        for (i = 1; i <= dim; i++) {
          ps++;
          p = v2;
          for (j = 1; j <= dim; j++) {
            *ps += (*(++p) * spm[j][i]);
            *ps %= prime;
          }
        }
      }
      else {
        p = v2;
        r = p + dim;
        if (mult)
          while (++p <= r) {
            *(++ps) += *p;
            x = cord[p - v2];
            *ps %= x;
          }
        else
          while (++p <= r) {
            *(++ps) += *p;
            *ps %= prime;
          }
      }
    }
    l = 0;
    ps = val;
    qs = val + dim;
    while (++ps <= qs)
      if (*ps != 0)
        l += 2;
    fprintf(op, "%4d  ", l);
    ps = val;
    while (++ps <= qs)
      if (*ps != 0)
        fprintf(op, "%4d%4d", (int)(ps - val), *ps);
    fprintf(op, "\n");
    if (mult == 0)
      fclose(ip);
  }
  fclose(ipx);
  unlink(outft);
  return (0);
}

int scan(int bno, short * la)
/* This is the computation of one word as described above */
{
  short cos, nl, min, canct, **cim;
  short c, *v1, *conco, *p, *q, *r, ***cco, *wp, *wr;
  v1 = wd1 + wv;
  cco = coeff[bno];
  cim = imcos[bno];
  cos = 1;
  wp = wd2;
  wr = wp + *la;
  nl = 0;
  while (++wp <= wr) {
    conco = cco[*wp][cos];
    cos = cim[*wp][cos];
    if (conco != 0) {
      q = conco + *conco + 1;
      r = q + *q;
      if (mult) {
        while (++q < r) {
          c = *(q++);
          v1[c] += *q;
          v1[c] %= cord[c];
        }
      }
      else {
        action(v1, conco);
        while (++q <= r) {
          c = *(q++);
          v1[c] += *q;
          v1[c] %= prime;
        }
      }
      min = (nl <= *conco) ? nl : *conco;
      canct = 0;
      p = wd1 + nl;
      q = conco + 1;
      while (canct < min && *p == invg[*q]) {
        canct++;
        p--;
        q++;
      }
      r = conco + *conco;
      while (q <= r)
        *(++p) = *(q++);
      nl += (*conco - 2 * canct);
    }
  }
  if (cos != 1) {
    fprintf(stderr,
            "One of the relations is not satisfied by the permutations.\n");
    return (-1);
  }
  if (bno > 1) {
    p = v1;
    q = wd2 + wv;
    r = p + dim;
    if (mult)
      while (++p <= r) {
        *p += *(++q);
        *p %= cord[p - v1];
        *q = 0;
      }
    else
      while (++p <= r) {
        *p += *(++q);
        *p %= prime;
        *q = 0;
      }
  }
  *la = nl;
  if (nl > wdl) {
    fprintf(stderr, "Word too long. Increase WDL.\n");
    return (-1);
  }
  return (0);
}

int action(short * a, short * b)
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
