#include "defs.h"

extern char inf0[], inf1[], inf4[], outf0[], outf1[], outf2[], outcopy[], act,
    ch1, norm;
extern short intexp, facexp, stage, mcl, prime, exp, nng, *rpf, *rpb, *eexpnt,
    *enexpnt, **pcb, dim, onng, **npcb, rel[], expnt[], nexpnt[], pinv[],
    wt[], d1[], d2[], *pcptr[], **powptr[], **comptr[], sd1[], sd2[], swt[],
    dpth[], **mat[], cp[], **intg, **cintg, cbno, ngens, maxm, matcl, **extno,
    **subno, chsdim, chpdim, exp1;
extern int   rsp, ptrsp, rspk, marg;
extern long  inf3offset, inf4offset;
short *      wf, *wc, **inteno, intcpd, inng, *nsd1, *nsd2, intfe, **npcb2;
extern FILE *ip, *ipm, *op;
FILE *       ipcopy, *opcopy;

void enter(short * g, int pow)
/* Enters a power of the word pointed to by g into rel for collection. */
{
  short *ps, *pf, *pc, i, sgn;
  if (pow < 0) {
    sgn = -1;
    pow = -pow;
    pf = g + 1;
    ps = g + *g + 1;
  }
  else {
    sgn = 1;
    ps = g - 1;
    pf = g + *g - 1;
  }
  for (i = 1; i <= pow; i++) {
    pc = ps;
    while (pc != pf) {
      pc += (2 * sgn);
      wc += 2;
      *wc = *pc;
      *(wc + 1) = *(pc + 1) * sgn;
    }
  }
}

void entvec(short * h, short * g, int pow)
/* Similar, but g is a vector. Used only for H^1. */
{
  short i, j, l, *p, *q;
  for (i = 1; i <= dim; i++)
    if ((j = g[i]) != 0) {
      j *= pow;
      j %= prime;
      if (j < 0)
        j += prime;
      wc += 2;
      *wc = exp + i;
      *(wc + 1) = j;
    }
  if (h != 0) {
    l = *h;
    p = h + l;
    while (++h < p) {
      q = nexpnt + (*h);
      *q += (pow * *(++h));
      *q %= prime;
      if (*q < 0)
        *q += prime;
    }
  }
}

void expand(short * p1, short * p2, int len)
/* Expand word p1 to vector p2 of length len. */
{
  short l, *p, *q;
  l = *p1;
  p = p1;
  q = p + l;
  zero(p2, p2 + len);
  while (++p < q) {
    p2[*p] = *(p + 1);
    p++;
  }
}

void compress(short * p1, short * p2, int len)
/* compress vector p2 of length len to word p1. */
{
  short i, n, *p;
  p = p1;
  for (i = 1; i <= len; i++)
    if ((n = p2[i]) != 0) {
      *(++p) = i;
      *(++p) = n;
    }
}

void nchg(void)
/* We have a new generator for H^2. This is put into echelon form with the
   existing generators, and stored.
*/
{
  short l, m, n, len, *p, *q, *v1, *r, *r1;
  for (m = 1; m <= onng; m++)
    if ((n = nexpnt[m]) != 0) {
      if (n < 0) {
        n += prime;
        nexpnt[m] = n;
      }
      if ((p = extno[m]) != 0) {
        expand(p, rpf, onng);
        p = rpf;
        n = prime - n;
        q = p + onng;
        v1 = nexpnt + 1;
        while (++p <= q) {
          *v1 += (n * *p);
          *v1 %= prime;
          v1++;
        }
      }
    }
  for (m = 1; m <= onng; m++)
    if ((n = nexpnt[m]) != 0) {
      n = pinv[n];
      p = nexpnt + m - 1;
      q = nexpnt + onng;
      while (++p <= q) {
        *p *= n;
        *p %= prime;
      }
      chpdim++;
      printf("chpdim=%2d,    bno=%2d\n", chpdim, m);
      for (l = 1; l <= onng; l++)
        if ((p = extno[l]) != 0) {
          r = rpb - onng;
          expand(p, r, onng);
          if ((n = r[m]) != 0) {
            n = prime - n;
            r1 = r;
            q = r + onng;
            v1 = nexpnt + 1;
            while (++r1 <= q) {
              *r1 += (n * *v1);
              *r1 %= prime;
              v1++;
            }
            len = 0;
            for (n = 1; n <= onng; n++)
              if (r[n] != 0)
                len += 2;
            if (len > *p) {
              p = rpf;
              extno[l] = p;
              rpf += (len + 1);
            }
            *p = len;
            compress(p, r, onng);
          }
        }
      extno[m] = rpf;
      len = 0;
      for (n = 1; n <= onng; n++)
        if (nexpnt[n] != 0)
          len += 2;
      *rpf = len;
      compress(rpf, nexpnt, onng);
      rpf += (len + 1);
      break;
    }
}

int spact(void)
/* Computes cohomology group H^i(P,M) or H^i(Q,M) */
{
  short i, j, k, l, m, ie, ct, cl, fg, wi, wj, *p, *q, *nrpf, *v1, *v2,
      **swop, homcl, c;
  char inp;
  inp = (ch1 && act) ? 0 : 1;
  /* inp as in calcfm in case ch1 */
  if (ingp(inp) == -1)
    return (-1);
  if (ch1 && act) {
    inf3offset = ftell(ip);
    fclose(ip);
  }
  ct = ngens;
  j = exp;
  for (i = 1; i <= dim; i++) {
    j++;
    wt[j] = swt[i];
    d1[j] = (sd1[i] > 0) ? sd1[i] + exp : 0;
    d2[j] = sd2[i];
  }
  printf("Computing matrices.\n");
  /* Matrices for all pcp gens of P or Q are now computed */
  if (maxm < 2 * facexp) {
    fprintf(stderr, "Not enough mat space. Increase MSP (of MV or MM).\n");
    return (-1);
  }
  for (i = facexp; i >= 1; i--)
    if (wt[i] == 1) {
      if (i > ct) {
        swop = mat[i];
        mat[i] = mat[ct];
        mat[ct] = swop;
        for (j = 1; j <= dim; j++)
          if (d2[exp + j] == ct)
            d2[exp + j] = i;
      }
      if (inv(mat[i], mat[facexp + i]) == -1)
        return (-1);
      ct--;
    }
  if (ct != 0) {
    fprintf(stderr, "No of pgens wrong.\n");
    return (-1);
  }
  for (i = 2; i <= facexp; i++)
    if (wt[i] > 1) {
      p = (d1[i] == d2[i]) ? *powptr[d1[i]] : *(comptr[d1[i]] + d2[i]);
      q = p + *p - 2;
      *cp = 0;
      while (--q > p) {
        k = *(q + 1);
        for (j = 1; j <= k; j++)
          cp[++(*cp)] = *q + facexp;
        q--;
      }
      if (d1[i] == d2[i])
        for (j = 1; j <= prime; j++)
          cp[++(*cp)] = d1[i];
      else {
        cp[++(*cp)] = d1[i] + facexp;
        cp[++(*cp)] = d2[i] + facexp;
        cp[++(*cp)] = d1[i];
        cp[++(*cp)] = d2[i];
      }
      prod(cp, mat[i]);
      if (inv(mat[i], mat[i + facexp]) == -1)
        return (-1);
    }
  if (act == 0) {
    op = fopen(outf2, "w");
    fprintf(op, "%3d %3d %3d\n", prime, dim, facexp);
    for (i = 1; i <= facexp; i++)
      printmat(mat[i]);
    fclose(op);
  }


  printf("Computing full pcp.\n");
  /* pcp of P or Q in dual action on module is computed from the matrices. */
  v1 = mat[facexp + 1][1];
  v2 = mat[facexp + 1][2];
  for (i = 1; i <= dim; i++)
    v1[i] = 0;
  ie = exp;
  for (i = 1; i <= dim; i++) {
    if (i > 1)
      v1[i - 1] = 0;
    v1[i] = 1;
    ie++;
    for (j = 1; j <= facexp; j++)
      if (comm(v1, v2, mat[j])) {
        *(comptr[ie] + j) = rpf + 1;
        nrpf = rpf + 2;
        l = 0;
        for (k = 1; k <= dim; k++)
          if ((m = v2[k]) != 0) {
            *(nrpf++) = k + exp;
            *(nrpf++) = m;
            l += 2;
          }
        *(rpf + 1) = l;
        m = *(nrpf - 2);
        if (d1[m] == ie && d2[m] == j)
          *rpf = m;
        else
          *rpf = 0;
        rpf = nrpf;
      }
  }
  if (ch1 == 0) {
    printf("Computing P-homs.\n");
    fflush(stdout);
    homcl = 0;
    /* Hom-P(FM,M) will now be computed  as dual of tensor product, using NQA
     */
    if (matcl + 1 >= mcl) {
      fprintf(stderr, "Class too big. Increase MCL.\n");
      return (-1);
    }
    for (cl = matcl + 1; cl > 1; cl--) {
      printf("cl=%d.\n", cl);
      for (fg = facexp + 1; dpth[fg] == 1 && fg <= exp; fg++)
        for (i = exp + 1; i <= exp + dim; i++)
          if (wt[i] + 1 == cl) {
            if (intgen(i, fg) == -1)
              return (-1);
            if ((k = wt[i] + wt[fg]) > homcl) {
              homcl = k;
              if (homcl + 1 >= mcl) {
                fprintf(stderr, "Class too big. Increase MCL.\n");
                return (-1);
              }
            }
          }
      while (fg <= exp) {
        for (i = exp + 1; i <= exp + dim; i++)
          if (wt[i] + dpth[fg] == cl) {
            if (dpth[d1[fg]] == dpth[fg] - 1) {
              if (assoc(i, d1[fg], d2[fg]))
                if (subrel(i, fg) == -1)
                  return (-1);
            }
            else if (intgen(i, fg) == -1)
              return (-1);
            if ((k = wt[i] + wt[fg]) > homcl) {
              homcl = k;
              if (homcl + 1 >= mcl) {
                fprintf(stderr, "Class too big. Increase MCL.\n");
                return (-1);
              }
            }
          }
        fg++;
      }
      for (i = 1; i <= facexp; i++) {
        wi = wt[i];
        for (j = facexp + 1; j <= exp; j++) {
          wj = dpth[j];
          if (wi + wj >= cl)
            break;
          for (k = exp + 1; k <= exp + dim; k++)
            if (wi + wj + wt[k] == cl)
              if (assoc(k, j, i)) {
                if ((l = prnrel()) == 0)
                  goto nextcl;
                if (l == -1)
                  return (-1);
              }
        }
      }
    nextcl:;
    }
    bgc();
    rsp = rpb - rel + 1;
    printf("Computing extensions. nng,homcl=%d,%d\n", nng, homcl);
    fflush(stdout);
    stage = 2;
    onng = nng;
    chpdim = 0;
    chsdim = 0;
    extno = pcptr + ptrsp - onng - 1;
    for (i = 1; i <= nng; i++)
      extno[i] = 0;
  }
  else {
    stage = 4;
    printf("Computing first cohomology group.\n");
    homcl = matcl;
    if (homcl + 1 >= mcl) {
      fprintf(stderr, "Class too big. Increase MCL.\n");
      return (-1);
    }
  }

  /* H^2 or H^1 will now be computed */
  for (cl = homcl + 1; cl >= 2; cl--) {
    printf("cl=%d\n", cl);
    if (cl <= matcl + 1) {
      for (i = 1; i <= facexp; i++)
        if (wt[i] == 1)
          for (j = exp + 1; j <= exp + dim; j++)
            if (wt[j] + 1 == cl)
              if (intgen(j, i) == -1)
                return (-1);
      for (i = 1; i <= facexp; i++) {
        wi = wt[i];
        if (wi > 1)
          for (j = exp + 1; j <= exp + dim; j++)
            if (wi + wt[j] == cl)
              if (assoc(j, d1[i], d2[i]))
                if (subrel(j, i) == -1)
                  return (-1);
      }
    }
    for (i = 1; i <= facexp; i++)
      for (j = i; j <= facexp; j++)
        for (k = exp + 1; k <= exp + dim; k++)
          if ((i == j && wt[i] + 1 + wt[k] == cl) ||
              (i != j && wt[i] + wt[j] + wt[k] == cl))
            if (assoc(k, j, i)) {
              if (ch1) {
                if ((l = prnrel()) == 0)
                  goto ncl;
                if (l == -1)
                  return (-1);
              }
              else {
                l = prnrel();
                if (l == -1)
                  return (-1);
                if (l > 1)
                  nchg();
              }
            }
  ncl:;
  }
  if (ch1)
    printf("Cohomology grp has dim %d\n", nng);
  /* We now no longer need the data for Q in the back of rel. */
  rpb = rel + rspk - 1;
  if (act) {
    if (ch1)
      onng = nng;
    else
      unlink(outf1);
  }
  else {
    if (stage > 2) {
      onng = nng;
      strcpy(outf1, outf0);
    }
    outgp();
    if (ch1) {
      ipcopy = fopen(outf1, "r");
      opcopy = fopen(outcopy, "w");
      while ((c = getc(ipcopy)) != -1)
        putc(c, opcopy);
      fclose(ipcopy);
      fclose(opcopy);
    }
  }
  return (0);
}

int intact(void)
/* In case act, appropriate quotient of H^i(P,M) is computed */
{
  short a, b, c, d, i, j, k, l, len, m, n, x, f1, f2, *p, *p1, *q, *q1, *r,
      *r1, *ia, *ib, *v1, **bim, **cbim, **imcg, **cimcg, *null;
  char pow, nz;
  strcpy(inf1, inf0);
  if (norm) {
    if (ch1)
      strcpy(inf1, outcopy);
    if (ingp(1) == -1)
      return (-1);
    if (ch1)
      strcpy(inf1, inf0);
    cbno = 1;
    bim = mat[cbno];
    onng = nng;
    for (i = 1; i <= dim; i++) {
      for (j = 1; j <= dim; j++)
        bim[i][j] = 0;
      bim[i][i] = 1;
    }
  }
  /* It is now time to read in the full pcp for the Frattini extension of P
     again, in order to compute the action on it. Before doing this, we must
     save all necessary information about Q. intg and cintg are already safely
     at the back of rel., and the cohomolgy gens of Q will also be copied
     there. The definitions of the gens of Q will be copied to sd1 and sd2,
     but a slightly messy method has been used to save long definitions.
  */
  p = d1;
  q = p + exp + dim + onng;
  v1 = sd1 + 1;
  while (++p <= q)
    *(v1++) = *p;
  nsd1 = sd1 + exp + dim;
  p = d2;
  q = p + exp + dim + onng;
  v1 = sd2 + 1;
  while (++p <= q)
    *(v1++) = *p;
  nsd2 = sd2 + exp + dim;
  if (ch1) {
    for (i = exp + 1; i <= exp + dim + nng; i++) {
      p = (sd1[i] == 0) ? 0 : *(comptr[sd1[i]] + sd2[i]);
      if (p != 0) {
        l = *p;
        if (i <= exp + dim)
          l -= 2;
        if (l > 0) {
          sd1[i] = -sd1[i];
          q = p + l;
          while (q > p) {
            *rpb = *q;
            q--;
            rpb--;
          }
          *rpb = l;
          cintg[i] = rpb;
          rpb--;
        }
      }
    }
  }
  else
    for (i = 1; i <= exp1; i++)
      if (sd1[i] != 0) {
        p = sd1[i] == sd2[i] ? *powptr[sd1[i]] : *(comptr[sd1[i]] + sd2[i]);
        l = *p;
        if (l > 2) {
          sd1[i] = -sd1[i];
          q = p + l - 2;
          while (q > p) {
            *rpb = *q;
            q--;
            rpb--;
          }
          *rpb = l - 2;
          cintg[i] = rpb;
          rpb--;
        }
      }
  inng = onng;
  intexp = exp;
  if (ch1 == 0) {
    p = wt;
    q = p + facexp;
    v1 = swt + 1;
    while (++p <= q)
      *(v1++) = *p;
    inteno = extno;
    intcpd = chpdim;
    intfe = facexp;
    if (norm == 0) {
      for (i = 1; i <= onng; i++)
        if ((p = inteno[i]) != 0) {
          q = p + *p;
          rpb -= (*p + 1);
          v1 = rpb;
          while (p <= q)
            *(++v1) = *(p++);
          inteno[i] = rpb + 1;
        }
      ptrsp -= onng;
    }
  }
  rsp = rpb - rel + 1;
  if (norm == 0 || ch1)
    if (ingp(1) == -1)
      return (-1);
  printf("Reading conj matrix.\n");
  /* matrix of dcrep is read and multiplied by action base change matrix */
  bim = mat[cbno];
  cbim = mat[cbno + 1];
  ip = fopen(inf4, "r");
  fseek(ip, inf4offset, 0);
  readmat(mat[cbno + 2]);
  inf4offset = ftell(ip);
  fclose(ip);
  *cp = 2;
  cp[1] = cbno;
  cp[2] = cbno + 2;
  prod(cp, cbim);
  if (ch1) {
    printf("Computing action.\n");
    wf = rpf;
    null = 0;
    for (i = 1; i <= dim; i++)
      *(++npcb2) = 0;
    if (npcb2 - pcptr >= ptrsp) {
      fprintf(stderr, "Out of ptrsp. Increase PTRSP.\n");
      return (-1);
    }
    for (i = intexp + 1; i <= intexp + dim + inng; i++) {
      wc = wf - 2;
      a = sd1[i];
      if (a == 0)
        continue;
      b = sd2[i];
      zero(nexpnt, enexpnt);
      if (a < 0) {
        p = cintg[i];
        q = p + *p - 1;
        while (q > p) {
          entvec(npcb[*q - intexp], bim[*q - intexp], -*(q + 1));
          q -= 2;
        }
      }
      ia = bim[abs(a) - intexp];
      ib = intg[b];
      entvec(null, ia, -1);
      enter(ib, -1);
      entvec(null, ia, 1);
      enter(ib, 1);
      zero(expnt, eexpnt);
      collect(wc, wf, 1);
      wc = wf - 2;
      if (a < 0) {
        p = cintg[i];
        q = p + *p - 1;
        while (q > p) {
          entvec(null, cbim[*q - intexp], -*(q + 1));
          q -= 2;
        }
      }
      ia = cbim[abs(a) - intexp];
      ib = cintg[b];
      entvec(null, ia, -1);
      enter(ib, -1);
      entvec(null, ia, 1);
      enter(ib, 1);
      zero(expnt, eexpnt);
      collect(wc, wf, -1);
      if (i <= intexp + dim) {
        l = 0;
        for (j = nng; j >= 1; j--)
          if ((k = nexpnt[j]) != 0) {
            if (k < 0)
              k += prime;
            *rpb = k;
            *(rpb - 1) = j;
            rpb -= 2;
            l += 2;
          }
        if (l != 0) {
          npcb[i - intexp] = rpb;
          *rpb = l;
          *(rpb - 1) = 0;
          rpb -= 2;
        }
      }
      else if ((l = prnrel()) == -1)
        return (l);
      else if (l == 0)
        return (2);
    }
    printf("End of this action. Present dimension=%d\n", nng);
    fflush(stdout);
    if (nng == 0)
      return (2);
  }
  else {
    printf("Computing Frattini embeddings.\n");
    /* Now intg and cintg are computed as strings in gens of P for the gens of
     * Q */
    for (i = 1; i <= exp1; i++)
      if (i > intfe || swt[i] > 1) {
        wf = rpf;
        a = sd1[i];
        b = sd2[i];
        wc = wf - 2;
        pow = abs(a) == b;
        /* a<0 means a long def */
        if (a < 0) {
          p = cintg[i];
          q = p + *p - 1;
          while (q > p) {
            enter(intg[*q], -*(q + 1));
            q -= 2;
          }
        }
        ia = intg[abs(a)];
        ib = intg[b];
        if (pow)
          enter(ia, prime);
        else {
          enter(ia, -1);
          enter(ib, -1);
          enter(ia, 1);
          enter(ib, 1);
        }
        zero(expnt, eexpnt);
        if (wc >= wf)
          collect(wc, wf, 1);
        l = 0;
        for (k = exp; k >= 1; k--)
          if ((x = expnt[k]) != 0) {
            *rpb = x;
            *(rpb - 1) = k;
            rpb -= 2;
            l += 2;
          }
        *rpb = l;
        intg[i] = rpb;
        rpb--;
        wc = wf - 2;
        if (a < 0) {
          a = -a;
          p = cintg[i];
          q = p + *p - 1;
          while (q > p) {
            enter(cintg[*q], -*(q + 1));
            q -= 2;
          }
        }
        ia = cintg[a];
        ib = cintg[b];
        if (pow)
          enter(ia, prime);
        else {
          enter(ia, -1);
          enter(ib, -1);
          enter(ia, 1);
          enter(ib, 1);
        }
        zero(expnt, eexpnt);
        if (wc >= wf)
          collect(wc, wf, 1);
        l = 0;
        for (k = exp; k >= 1; k--)
          if ((x = expnt[k]) != 0) {
            *rpb = x;
            *(rpb - 1) = k;
            rpb -= 2;
            l += 2;
          }
        *rpb = l;
        cintg[i] = rpb;
        rpb--;
        if (rpb - rpf < marg) {
          fprintf(stderr, "Out of space. Increase RSP.\n");
          return (-1);
        }
      }
    printf("Computing images of cohomology gens.\n");
    /* Now images of cohomology gens of Q and its conjugte are computed as
       strings, and stored in imcg and cimcg.
    */

    imcg = pcb;
    cimcg = imcg + inng;
    if (pcb + 2 * (inng + nng) - pcptr >= ptrsp) {
      fprintf(stderr, "Out of ptr space. Increase PTRSP.\n");
      return (-1);
    }
    for (i = 1; i <= inng; i++) {
      a = nsd1[i] - intexp;
      b = nsd2[i];
      zero(rpf, rpf + nng);
      p = bim[a];
      p1 = p + dim + 1;
      q = intg[b];
      q1 = q + *q;
      while (--p1 > p)
        if ((f1 = *p1) != 0) {
          c = p1 - p + exp;
          q = intg[b];
          while (++q < q1) {
            d = *(q++);
            f2 = *q;
            if ((r = *(comptr[c] + d)) != 0) {
              r1 = r + *r;
              while (++r < r1) {
                rpf[*r] += (*(r + 1) * f1 * f2);
                r++;
              }
            }
          }
        }
      p = rpf;
      p1 = p + nng;
      while (++p <= p1)
        *p %= prime;
      len = 0;
      for (j = 1; j <= nng; j++)
        if (rpf[j] != 0)
          len += 2;
      rpb -= len;
      imcg[i] = rpb;
      *rpb = len;
      compress(rpb, rpf, nng);
      rpb--;
      zero(rpf, rpf + nng);
      p = cbim[a];
      p1 = p + dim + 1;
      q = cintg[b];
      q1 = q + *q;
      while (--p1 > p)
        if ((f1 = *p1) != 0) {
          c = p1 - p + exp;
          q = cintg[b];
          while (++q < q1) {
            d = *(q++);
            f2 = *q;
            if ((r = *(comptr[c] + d)) != 0) {
              r1 = r + *r;
              while (++r < r1) {
                rpf[*r] += (*(r + 1) * f1 * f2);
                r++;
              }
            }
          }
        }
      p = rpf;
      p1 = p + nng;
      while (++p <= p1)
        *p %= prime;
      len = 0;
      for (j = 1; j <= nng; j++)
        if (rpf[j] != 0)
          len += 2;
      rpb -= len;
      cimcg[i] = rpb;
      *rpb = len;
      compress(rpb, rpf, nng);
      rpb--;
      if (rpb - rpf - marg > rsp) {
        fprintf(stderr, "Out of space. Increase RSP.\n");
        return (-1);
      }
    }
    printf("Computing Sub-cohomology generators.\n");
    /* Finally we are ready to compute gens of the subgroup of H^2(P,M) to be
       factored out. These will be put into echelon form.
    */
    for (i = 1; i <= inng; i++)
      if ((p = inteno[i]) != 0) {
        zero(nexpnt, enexpnt);
        p1 = p + *p;
        while (--p1 > p) {
          a = *p1;
          f1 = *(p1 + 1);
          p1--;
          q = rpf;
          q1 = q + nng;
          r = q1;
          r1 = nexpnt;
          expand(imcg[a], q, nng);
          expand(cimcg[a], r, nng);
          while (++q <= q1) {
            r++;
            r1++;
            *r1 += (f1 * *q);
            *r1 -= (f1 * *r);
          }
        }
        r = nexpnt;
        nz = 0;
        r1 = r + nng;
        while (++r <= r1) {
          *r %= prime;
          if (*r != 0)
            nz = 1;
        }
        if (nz) {
          for (m = 1; m <= nng; m++)
            if ((n = nexpnt[m]) != 0) {
              if (n < 0) {
                n += prime;
                nexpnt[m] = n;
              }
              if ((q = subno[m]) != 0) {
                n = prime - n;
                r = q + *q;
                while (++q < r) {
                  q1 = nexpnt + *q;
                  q++;
                  *q1 += (n * *q);
                  *q1 %= prime;
                }
              }
            }
          for (m = 1; m <= nng; m++)
            if ((n = nexpnt[m]) != 0) {
              n = pinv[n];
              q = nexpnt + m - 1;
              r = nexpnt + nng;
              while (++q <= r) {
                *q *= n;
                *q %= prime;
              }
              chsdim++;
              printf("chsdim=%2d,    bno=%2d\n", chsdim, m);
              for (l = 1; l <= nng; l++)
                if ((p = subno[l]) != 0) {
                  r = rpb - nng;
                  expand(p, r, nng);
                  if ((n = r[m]) != 0) {
                    n = prime - n;
                    q = r + nng;
                    q1 = nexpnt + 1;
                    r1 = r;
                    while (++r1 <= q) {
                      *r1 += (n * *q1);
                      *r1 %= prime;
                      q1++;
                    }
                    len = 0;
                    for (n = 1; n <= nng; n++)
                      if (r[n] != 0)
                        len += 2;
                    if (len > *p) {
                      p = rpf;
                      subno[l] = p;
                      rpf += (len + 1);
                    }
                    *p = len;
                    compress(p, r, nng);
                  }
                }
              subno[m] = rpf;
              len = 0;
              for (n = 1; n <= nng; n++)
                if (nexpnt[n] != 0)
                  len += 2;
              *rpf = len;
              compress(rpf, nexpnt, nng);
              rpf += (len + 1);
              if (rpb - rpf < marg) {
                fprintf(stderr, "Out of space. Increase RSP.\n");
                return (-1);
              }
              break;
            }
        }
      }
    printf("End of this action. chpdim,chsdim=%d,%d\n\n", chpdim, chsdim);
    fflush(stdout);
    if (chpdim == chsdim)
      return (2);
  }
  strcpy(outf1, inf0);
  onng = nng;
  outgp();
  return (0);
}
