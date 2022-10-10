#include "defs.h"

extern char  inf1[], outf[], outfm[], gap;
extern short intexp, mexp, mng, wksp, prime, exp, nng, class, *rpf, *rpb,
    *eexpnt, *enexpnt, **pcb, mnng, mord, rel[], expnt[], nexpnt[], cord[],
    wt[], d1[], d2[], *pcptr[], **powptr[], **comptr[], *sspc[], *sspf[],
    sgen[], sex[], spgen[], spex[], spugen[], *tlintg[];
extern int ptrsp, rsp;
short      fac;

int ingp(void)
/* The PCP for P (output from pcrun) is read in, and all of the associated
   pointers are set. If -a is true, then the input comes from a preceding run
   of nqrun, and tails are also read in.
*/
{
  short  i, j, k, l, m, x, y, no, *orpf, *orpb, **pcp;
  char   tails;
  FILE * ip;
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd%hd%hd", &prime, &exp, &nng, &no, &class, &m);
  if (exp > mexp) {
    fprintf(stderr, "exp too big. Increase MEXP.\n");
    return (-1);
  }
  if (nng > mng) {
    fprintf(stderr, "nng too big. Increase MNG.\n");
    return (-1);
  }
  if (m != 1) {
    fprintf(stderr, "Wrong version of nq.\n");
    return (-1);
  }
  tails = exp == no;
  if (tails == 0)
    nng = 0;
  for (i = 1; i <= exp; i++)
    fscanf(ip, "%hd", wt + i);
  for (i = 1; i <= exp + nng; i++)
    fscanf(ip, "%hd", d1 + i);
  for (i = 1; i <= exp + nng; i++)
    fscanf(ip, "%hd", d2 + i);
  for (i = 1; i <= nng; i++)
    fscanf(ip, "%hd", cord + i);
  rpf = rel;
  rpb = rel + rsp - 1;
  pcp = pcptr;
  for (i = 2; i <= no; i++) {
    comptr[i] = pcp - 2;
    for (j = 1; j < i; j++) {
      fscanf(ip, "%hd%hd", &k, &l);
      if (l == 0)
        *(pcp++) = 0;
      else {
        *(pcp++) = rpf + 1;
        *(rpf++) = k;
        *rpf = l;
        orpf = rpf;
        rpf += (l + 1);
        while ((++orpf) < rpf)
          fscanf(ip, "%hd", orpf);
      }
      if (tails) {
        fscanf(ip, "%hd%hd", &k, &l);
        if (l == 0)
          *(pcp++) = 0;
        else {
          orpb = rpb - nng;
          *(pcp++) = orpb;
          *orpb = k;
          while (rpb > orpb)
            *(rpb--) = 0;
          l /= 2;
          for (k = 1; k <= l; k++) {
            fscanf(ip, "%hd%hd", &x, &y);
            rpb[x] = y;
          }
          rpb--;
        }
      }
      else
        *(pcp++) = 0;
    }
  }
  for (i = no + 1; i <= exp; i++) {
    comptr[i] = pcp - 2;
    for (j = 1; j < i; j++) {
      *(pcp++) = 0;
      *(pcp++) = 0;
    }
  }
  for (i = 1; i <= no; i++) {
    powptr[i] = pcp;
    fscanf(ip, "%hd%hd", &k, &l);
    if (l == 0)
      *(pcp++) = 0;
    else {
      *(pcp++) = rpf + 1;
      *(rpf++) = k;
      *rpf = l;
      orpf = rpf;
      rpf += (l + 1);
      while ((++orpf) < rpf)
        fscanf(ip, "%hd", orpf);
    }
    if (tails) {
      fscanf(ip, "%hd%hd", &k, &l);
      if (l == 0)
        *(pcp++) = 0;
      else {
        orpb = rpb - nng;
        *(pcp++) = orpb;
        *orpb = k;
        while (rpb > orpb)
          *(rpb--) = 0;
        l /= 2;
        for (k = 1; k <= l; k++) {
          fscanf(ip, "%hd%hd", &x, &y);
          rpb[x] = y;
        }
        rpb--;
      }
    }
    else
      *(pcp++) = 0;
  }
  for (i = no + 1; i <= exp; i++) {
    powptr[i] = pcp;
    *(pcp++) = 0;
    *(pcp++) = 0;
  }
  pcb = pcp - 1;
  if (pcb - pcptr >= ptrsp) {
    fprintf(stderr, "Out of ptr space. Increase PTRSP.\n");
    return (-1);
  }
  if (rpb - rpf < wksp) {
    fprintf(stderr, "Out of space. Increase RSP.\n");
    return (-1);
  }
  fclose(ip);
  return (0);
}

void outgp(void)
/* The PCP is output, together with tails */
{
  short  i, k, l, **pcp, *b, *e, *c;
  FILE * op;
  op = fopen(outf, "w");
  fprintf(op, "%4d%4d%4d%4d%4d%4d\n", prime, exp, nng, exp, class, 1);
  for (i = 1; i <= exp; i++)
    fprintf(op, "%4d", wt[i]);
  fprintf(op, "\n");
  for (i = 1; i <= exp + nng; i++)
    fprintf(op, "%4d", d1[i]);
  fprintf(op, "\n");
  for (i = 1; i <= exp + nng; i++)
    fprintf(op, "%4d", d2[i]);
  fprintf(op, "\n");
  if (nng > 0) {
    for (i = 1; i <= nng; i++)
      fprintf(op, "%4d", cord[i]);
    fprintf(op, "\n");
  }
  pcp = pcptr;
  while (pcp < pcb) {
    if (*pcp == 0) {
      fprintf(op, "%3d  %3d", 0, 0);
      l = 0;
    }
    else {
      l = **pcp;
      fprintf(op, "%3d  %3d", *(*pcp - 1), l);
      b = *pcp;
      e = b + l;
      while (++b <= e)
        fprintf(op, "%4d", *b);
    }
    if (l < 10) {
      l = 10 - l;
      for (k = 1; k <= l; k++)
        fprintf(op, "    ");
    }
    pcp++;
    if (*pcp == 0)
      fprintf(op, "%3d  %3d", 0, 0);
    else {
      fprintf(op, "%3d  ", **pcp);
      b = *pcp;
      e = b + nng;
      c = cord + 1;
      l = 0;
      while (++b <= e) {
        while (*b < 0)
          (*b) += *c;
        if (*b != 0)
          l += 2;
        c++;
      }
      fprintf(op, "%3d", l);
      b = *pcp;
      while (++b <= e)
        if (*b != 0)
          fprintf(op, "%4d%4d", (int)(b - (*pcp)), *b);
    }
    fprintf(op, "\n");
    pcp++;
  }
  fclose(op);
}

void zero(short * p1, short * p2)
{
  while ((++p1) <= p2)
    *p1 = 0;
}

void setnr(short * p)
/* This is called only by collect, which follows */
{
  short *p1, *p2, *c;
  p1 = nexpnt;
  p2 = p1 + nng;
  c = cord + 1;
  while ((++p1) <= p2) {
    p++;
    *p1 += fac * (*p);
    *p1 %= *c;
    c++;
  }
}

int collect(short * spc, short * spf, int sgn)
/* The basic collecting algorithm. The word to be collected is stored as a
   gen-pow string. spc points to the last gen, and spf to the first, in this
   string. The collected word is put into expnt, as a vector, and its tail
   is added to nexpnt (or subtracted if sgn= -1). The algorithm is similar to
   that used in the Canberra NQA, but negative exponents are allowed.
*/
{
  short gen, ex, pgen, pex, pugen, stkp, exsgn, *e, *p1, *p2, *p3, *p4, **dp;
  stkp = 0;
  pgen = 0;
  pex = 0;
  pugen = 0;
  exsgn = 1;

recurse:
  gen = *spc;
  ex = exsgn * (*(spc + 1));
  stkp++;
  sspc[stkp] = spc;
  sspf[stkp] = spf;
  sgen[stkp] = gen;
  sex[stkp] = ex;
  spex[stkp] = pex;
  spugen[stkp] = pugen;
  spgen[stkp] = pgen;

loop:
  e = expnt + pgen;
  if (pgen == gen) {
    if (pugen == gen) {
      *e += ex;
      ex = 1;
      sex[stkp] = 1;
    }
    else
      (*e)++;
    if (*e < 0) {
      *e += prime;
      dp = powptr[gen];
      p1 = *dp;
      p2 = *(dp + 1);
      if (p2 != 0) {
        fac = -sgn;
        setnr(p2);
      }
      if (p1 != 0) {
        exsgn = -1;
        pex = -1;
        spc = p1 + 1;
        spf = spc + *p1 - 2;
        goto recurse;
      }
    }
    else if (*e >= prime) {
      *e -= prime;
      dp = powptr[gen];
      p1 = *dp;
      p2 = *(dp + 1);
      if (p2 != 0) {
        fac = sgn;
        setnr(p2);
      }
      if (p1 != 0) {
        exsgn = 1;
        pex = -1;
        pugen = gen;
        spf = p1 + 1;
        spc = spf + *p1 - 2;
        goto recurse;
      }
    }
  }
  else if (pex <= 0) {
    if (pugen == pgen)
      pugen++;
    pgen++;
    pex = *(e + 1);
    goto loop;
  }
  else {
    dp = comptr[gen] + 2 * pgen;
    p1 = *dp;
    p2 = *(dp + 1);
    if (p1 == 0) {
      if (p2 != 0) {
        fac = (pugen == pgen) ? ex * pex * sgn : pex * sgn;
        setnr(p2);
      }
      pex = 0;
      goto loop;
    }
    else {
      if (pugen == pgen) {
        spugen[stkp] = -pgen;
        if (ex < 0) {
          sex[stkp] = ex + prime;
          dp = powptr[gen];
          p3 = *dp;
          p4 = *(dp + 1);
          if (p4 != 0) {
            fac = -sgn;
            setnr(p4);
          }
          if (p3 != 0) {
            exsgn = -1;
            spc = p3 + 1;
            spf = spc + *p3 - 2;
            goto recurse;
          }
        }
      }
      else
        pugen = pgen;
    }
    spf = p1 + 1;
    spc = spf + *p1 - 2;
    exsgn = 1;
    if (p2 != 0) {
      fac = sgn;
      setnr(p2);
    }
    pex--;
    goto recurse;
  }

powdone:
  pgen = spgen[stkp];
  spc = sspc[stkp];
  spf = sspf[stkp];
  pex = spex[stkp];
  pugen = abs(spugen[stkp]);
  ex--;
  if (ex == 0) {
    if (spc == spf) {
      stkp--;
      if (stkp == 0)
        return (0);
      gen = sgen[stkp];
      ex = sex[stkp];
      if (pex == -1)
        goto powdone;
      spc = sspc[stkp];
      spf = sspf[stkp];
      pugen = spugen[stkp];
      goto loop;
    }
    else {
      if (spc < spf) {
        spc += 2;
        gen = *spc;
        ex = -*(spc + 1);
      }
      else {
        spc -= 2;
        gen = *spc;
        ex = *(spc + 1);
      }
      sgen[stkp] = gen;
      sspc[stkp] = spc;
      sex[stkp] = ex;
      pugen = pgen;
      goto loop;
    }
  }
  else {
    if (pugen > pgen) {
      pgen = pugen;
      pex = expnt[pgen];
    }
    sex[stkp] = ex;
    goto loop;
  }
}

int intgen(int i, int j)
/* A new generator is introduced into the tail of the relation [i,j], or
   i^p if i=j.
*/
{
  short **dp, *p, *nrpb, sum;
  if (i == j)
    dp = powptr[i];
  else
    dp = comptr[i] + 2 * j;
  if (*dp != 0 && *(*dp - 1) != 0) {
    printf("Reln [%d,%d] is already a defn.\n", i, j);
    return (0);
  }
  p = *(dp + 1);
  if (p != 0) {
    fprintf(stderr, "Reln [%d,%d] already has a tail.\n", i, j);
    return (-1);
  }
  nng++;
  sum = exp + nng;
  if (sum > mng) {
    fprintf(stderr, "Too many newgens. Increase MNG.\n");
    return (-1);
  }
  d1[sum] = i;
  d2[sum] = j;
  cord[nng] = mord;
  enexpnt++;
  nrpb = rpb - mnng;
  if (nrpb - rpf < wksp) {
    fprintf(stderr, "Out of space. Increase RSP.\n");
    return (-1);
  }
  zero(nrpb, rpb);
  *(dp + 1) = nrpb;
  nrpb[nng] = 1;
  *nrpb = sum;
  rpb = nrpb - 1;
  return (0);
}

int subrel(int i, int j)
/* The element currently in nexpnt, is introduced as the tail of the relation
   [i,j]
*/
{
  short **dp, *p, *nrpb;
  if (i == j)
    dp = powptr[i];
  else
    dp = comptr[i] + 2 * j;
  if (*dp != 0 && *(*dp - 1) != 0) {
    printf("Reln [%d,%d] is already a defn.\n", i, j);
    return (0);
  }
  p = *(dp + 1);
  if (p != 0) {
    fprintf(stderr, "Reln [%d,%d] already has a tail.\n", i, j);
    return (-1);
  }
  nrpb = rpb - mnng;
  if (nrpb - rpf < wksp) {
    fprintf(stderr, "Out of space. Increase RSP.\n");
    return (-1);
  }
  zero(nrpb, rpb);
  *(dp + 1) = nrpb;
  *nrpb = 0;
  rpb = nrpb - 1;
  for (i = 1; i <= nng; i++)
    nrpb[i] = nexpnt[i];
  return (0);
}

int assoc(int g1, int g2, int g3)
/* The difference between the tails of the elements (g1g2)g3 and g1(g2g3) is
   computed into nexpnt. Thisdifference will of course have to be made trivial
   later. If g1=g2, we use (g1^p)g2 and g1(g1^(p-1)g2), and similarly if
   g2=g3.
*/
{
  char  eq12, eq23, prnt, triv;
  short i, e, *p;
  if (g3 < 0) {
    prnt = 1;
    g3 = -g3;
  }
  else
    prnt = 0;
  eq12 = g1 == g2;
  eq23 = g2 == g3 && g1 != g2;
  zero(nexpnt, enexpnt);
  zero(expnt, eexpnt);
  p = rpf;
  *(p++) = g1;
  *(p++) = eq12 ? prime - 1 : 1;
  *p = g2;
  *(p + 1) = 1;
  collect(p, rpf, 1);
  p = rpf;
  for (i = 1; i <= exp; i++) {
    e = expnt[i];
    if (e != 0) {
      *(p++) = i;
      *(p++) = e;
    }
  }
  if (p != rpf) {
    zero(expnt, eexpnt);
    expnt[g3] = eq23 ? prime - 1 : 1;
    collect(p - 2, rpf, 1);
  }
  zero(expnt, eexpnt);
  if (eq12 && g2 == g3)
    expnt[g2] = 2;
  else {
    p = rpf;
    *(p++) = g2;
    *(p++) = 1;
    *p = g3;
    *(p + 1) = eq23 ? prime - 1 : 1;
    collect(p, rpf, -1);
  }
  p = rpf;
  *p = g1;
  *(p + 1) = eq12 ? prime - 1 : 1;
  collect(p, rpf, -1);
  if (prnt) {
    p = nexpnt;
    while (++p <= enexpnt)
      printf("%3d", *p);
    printf("\n");
  }
  triv = 1;
  for (i = 1; i <= nng; i++)
    if (triv && nexpnt[i] != 0)
      triv = 0;
  if (triv)
    return (0);
  else
    return (1);
}

int prnrel(int corrtl)
/* The element stored in nexpnt is made trivial. This involves eliminating or
   reducing the order of one or more new generators. corrtl is true when -a is
   set, and images of generators in the Sylow intersection have their own
   tails, which require updating. When orders involve larger powers of primes,
   sometimes more than one run is necessary to complete the reductions. Hence
   label "restart".
*/
{
  short i, j, k, len, hp, ct, fac, x, y, a, b, c, hpi, gno, *p, **dp,
      term[30];
  char elim, sub, rep, triv;
  a = 0;
  b = 0;
  gno = 0;
restart:
  hp = exp / 2;
  len = 0;
  rep = 0;
  for (i = nng; i >= 1; i--) {
    hpi = 0;
    x = nexpnt[i];
    if (x != 0) {
      len++;
      while (x % prime == 0) {
        hpi++;
        x /= prime;
      }
      if (hpi < hp || (hpi == hp && cord[i] < b)) {
        gno = i;
        hp = hpi;
        a = x;
        b = cord[i];
      }
    }
  }
  if (len == 0)
    return (0);
  elim = hp == 0;
  sub = len > 1;
  x = 1;
  for (i = 1; i <= hp; i++)
    x *= prime;
  if (sub) {
    ct = 0;
    c = cord[gno] / x;
    if (a < 0)
      a += cord[gno];
    while (a != 1) {
      ct++;
      term[ct] = c / a;
      b = c % a;
      c = a;
      a = b;
    }
    if (ct == 0)
      fac = 1;
    else {
      term[ct + 1] = 1;
      term[ct + 2] = 0;
      p = term + ct + 1;
      while (--p > term)
        *p = -(*p) * *(p + 1) + *(p + 2);
      fac = term[1];
    }
    p = nexpnt;
    while (++p <= enexpnt) {
      *p *= fac;
      *p /= x;
      *p %= cord[p - nexpnt];
    }
    dp = pcptr;
    while (++dp <= pcb) {
      p = *dp;
      if (p != 0) {
        y = p[gno];
        if (y != 0) {
          for (j = 1; j <= nng; j++)
            if (j == gno)
              continue;
            else if ((k = nexpnt[j]) != 0) {
              p[j] -= (y * k);
              p[j] %= cord[j];
            }
        }
      }
      dp++;
    }
    if (corrtl) {
      for (i = 3; i <= intexp; i++)
        if ((p = tlintg[i]) != 0) {
          y = p[gno];
          if (y != 0) {
            for (j = 1; j <= nng; j++)
              if (j == gno)
                continue;
              else if ((k = nexpnt[j]) != 0) {
                p[j] -= (y * k);
                p[j] %= cord[j];
              }
          }
        }
    }
    a = cord[gno];
    p = nexpnt;
    while (++p <= enexpnt) {
      *p *= a;
      *p %= cord[p - nexpnt];
      if (*p != 0)
        rep = 1;
    }
  }
  if (elim) {
    printf("New gen no %d is eliminated.\n", gno);
    dp = pcptr;
    while (++dp <= pcb) {
      p = *dp;
      if (p != 0) {
        if (*p > (exp + gno))
          (*p)--;
        else if (*p == (exp + gno))
          *p = 0;
        triv = 1;
        for (j = 1; j < gno; j++)
          if (p[j] != 0)
            triv = 0;
        for (j = gno; j < nng; j++) {
          p[j] = p[j + 1];
          if (triv && p[j] != 0)
            triv = 0;
        }
        p[nng] = 0;
        if (triv)
          *dp = 0;
      }
      dp++;
    }
    for (j = gno; j < nng; j++) {
      cord[j] = cord[j + 1];
      d1[exp + j] = d1[exp + j + 1];
      d2[exp + j] = d2[exp + j + 1];
      nexpnt[j] = nexpnt[j + 1];
    }
    if (corrtl)
      for (i = 3; i <= intexp; i++)
        if ((p = tlintg[i]) != 0)
          for (j = gno; j < nng; j++)
            p[j] = p[j + 1];
    nng--;
    mnng--;
    enexpnt--;
    if (nng == 0) {
      if (gap) {
        FILE * op = fopen(outfm, "w");
        fprintf(op, "COHOMOLO.Multiplier:=[];\n");
        fclose(op);
        printf("All new generators eliminated. Multiplier is trivial.\n");
      }
      else
        fprintf(stderr,
                "All new generators eliminated. Multiplier is trivial.\n");
      return (-1);
    }
  }
  else {
    printf("New order of new gen no %d is %d.\n", gno, x);
    cord[gno] = x;
    dp = pcptr;
    while (++dp <= pcb) {
      p = *dp;
      if (p != 0)
        p[gno] %= x;
      dp++;
    }
    if (corrtl)
      for (i = 3; i <= intexp; i++)
        if ((p = tlintg[i]) != 0)
          p[gno] %= x;
  }
  if (rep)
    goto restart;
  return (1);
}
