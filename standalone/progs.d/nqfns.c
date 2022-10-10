#include "defs.h"

extern char  inf[], inf1[], outf1[];
extern short facexp, tails, stage, depth, no, mng, mexp, prime, exp, nng,
    class, dim, onng, *rpf, *rpb, *eexpnt, *enexpnt, **pcb, **opcb, **npcb,
    **npcb2, *nd1, *nd2, **extno, **subno, chpdim, chsdim, rel[], expnt[],
    nexpnt[], prvec[], pinv[], wt[], d1[], d2[], *pcptr[], **powptr[],
    **comptr[], *sspc[], *sspf[], sgen[], sex[], spgen[], spex[], spugen[],
    dpth[];
extern int   rsp, wsp, ptrsp, marg;
short        fac;
extern FILE *ip, *op;

int ingp(int inp)
/* Input a group with pcp. Normally inp=1, and file inf1 is opened inside
   procedure. inp=0 when file is already open - i.e. when act=1, and pcp is
   the output from scrun.
*/
{
  int    i, j, k, l, sum, jump;
  short *orpf, **pcp, **pcpj, m;
  if (inp) {
    ip = fopen(inf1, "r");
    if (ip == 0) {
      fprintf(stderr, "Cannot open %s\n", inf1);
      return (-1);
    }
    fscanf(ip, "%hd%hd%hd%hd%hd%hd", &prime, &exp, &facexp, &no, &class, &m);
    if (exp >= mexp) {
      fprintf(stderr, "exp too big. Increase MEXP.\n");
      return (-1);
    }
    if (m != 0) {
      fprintf(stderr, "Wrong version of nq.\n");
      return (-1);
    }
  }
  nng = 0;
  sum = (stage == 2 || stage == 4) ? exp + dim : exp;
  if (inp)
    for (i = 1; i <= sum; i++)
      fscanf(ip, "%hd", wt + i);
  if (stage < 2) {
    if (exp > facexp)
      for (i = 1; i <= exp; i++)
        fscanf(ip, "%hd", dpth + i);
    else
      for (i = 1; i <= exp; i++)
        dpth[i] = 0;
    depth = dpth[exp];
  }
  for (i = 1; i <= sum; i++)
    fscanf(ip, "%hd", d1 + i);
  for (i = 1; i <= sum; i++)
    fscanf(ip, "%hd", d2 + i);
  /* Now we start to read the pcp into rel, setting the commutator and power
     pointers comptr and powptr, as we go.
     Defs and powers [d[j],d[i]] and d[j]^p are only in the file for j<=no, so
     for the remaining ones we zero them.
     pcp is the current place of the pointer in array of pointers pcptr.
     opcb,pcb,npcb and npcb2 are externals used to record the values of pcp
     at various key places, for use in outgp and prnrel.
     Roughly, opcb=value of pcp at end of PCP of P itself. (=pcb when
     stage=0). pcb=value at end of definition of commutators in [d(M),P]
     (d(M)=dual of M).
     npcb=value at end. npcb2 is only used for ch1.
  */
  rpf = rel;
  rpb = rel + rsp - 1;
  pcp = pcptr;
  for (i = 2; i <= no; i++) {
    comptr[i] = pcp - 1 - tails;
    m = (i > facexp) ? facexp : i - 1;
    for (j = 1; j <= m; j++) {
      fscanf(ip, "%d%d", &k, &l);
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
      if (tails)
        *(pcp++) = 0;
    }
  }
  if (tails || no < facexp)
    for (i = no + 1; i <= exp; i++) {
      comptr[i] = pcp - 1 - tails;
      m = (i > facexp) ? facexp : i - 1;
      for (j = 1; j <= m; j++) {
        *(pcp++) = 0;
        if (tails)
          *(pcp++) = 0;
      }
    }
  m = (no >= facexp) ? facexp : no;
  for (i = 1; i <= m; i++) {
    powptr[i] = pcp;
    fscanf(ip, "%d%d", &k, &l);
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
    if (tails)
      *(pcp++) = 0;
  }
  for (i = no + 1; i <= facexp; i++) {
    powptr[i] = pcp;
    *(pcp++) = 0;
    if (tails)
      *(pcp++) = 0;
  }
  if (stage == 3)
    no = exp;
  if (stage == 0)
    pcb = pcp - 1;
  /* End of input of PCP of P */
  else {
    opcb = pcp - 1;
    nd1 = d1 + exp + dim;
    nd2 = d2 + exp + dim;
    if (stage == 2 || stage == 4) {
      fscanf(ip, "%hd", &nng);
      for (i = 1; i <= nng; i++)
        fscanf(ip, "%hd", nd1 + i);
      for (i = 1; i <= nng; i++)
        fscanf(ip, "%hd", nd2 + i);
      if (stage == 4) {
        jump = exp * dim;
        pcpj = pcp + jump;
      }
      else {
        jump = 0;
        pcpj = NULL;
      }
    }
    else {
      jump = 0;
      pcpj = NULL;
    }
    for (i = exp + 1; i <= exp + dim; i++) {
      comptr[i] = pcp - 1;
      if (stage == 4)
        powptr[i] = pcpj - 1;
      for (j = 1; j <= exp; j++)
        if (stage == 3 || stage == 1)
          *(pcp++) = 0;
        else {
          fscanf(ip, "%d%d", &k, &l);
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
          if (stage == 4) {
            fscanf(ip, "%d%d", &k, &l);
            if (l == 0)
              *(pcpj++) = 0;
            else {
              orpf = rpb - l;
              *(pcpj++) = orpf;
              *(orpf - 1) = k;
              *orpf = l;
              while ((++orpf) <= rpb)
                fscanf(ip, "%hd", orpf);
              rpb -= (l + 2);
            }
          }
        }
    }
    pcb = pcp - 1;
    if (stage != 4 && stage != 2)
      for (i = exp + 1; i <= exp + dim; i++) {
        powptr[i] = pcp - 1;
        for (j = 1; j <= facexp; j++)
          *(pcp++) = 0;
      }
    else if (stage == 4)
      pcp += jump;
    npcb = pcp - 1;
    npcb2 = npcb;
  }
  wsp = 0;
  eexpnt = expnt + exp;
  if (stage)
    eexpnt += dim;
  enexpnt = nexpnt + nng;
  if (stage == 2) {
    extno = pcptr + ptrsp - 1 - nng;
    subno = extno - nng;
    for (i = 1; i <= nng; i++) {
      extno[i] = 0;
      subno[i] = 0;
    }
    fscanf(ip, "%hd", &chpdim);
    for (i = 1; i <= chpdim; i++) {
      fscanf(ip, "%d", &j);
      extno[j] = rpf;
      fscanf(ip, "%hd", rpf);
      l = *rpf;
      rpf++;
      for (j = 1; j <= l; j++) {
        fscanf(ip, "%hd", rpf);
        rpf++;
      }
    }
    fscanf(ip, "%hd", &chsdim);
    for (i = 1; i <= chsdim; i++) {
      fscanf(ip, "%d", &j);
      subno[j] = rpf;
      fscanf(ip, "%hd", rpf);
      l = *rpf;
      rpf++;
      for (j = 1; j <= l; j++) {
        fscanf(ip, "%hd", rpf);
        rpf++;
      }
    }
  }
  if (pcp - pcptr > ptrsp) {
    fprintf(stderr, "Not enough pointer space. Increase PTRSP.\n");
    return (-1);
  }
  if (rpf - rel > rsp) {
    fprintf(stderr, "Not enough space. Increase RSP.\n");
    return (-1);
  }
  if (inp)
    fclose(ip);
  return (0);
}

int outgp(void)
/* A group is output to outf1 */
{
  int     i, k, l, sum, nexp, jump;
  short **pcp, *b, *c, *e, *f, **pcpj;
  char    jtl;
  op = fopen(outf1, "w");
  sum = (stage) ? exp + dim : (tails) ? exp + nng : exp;
  if (tails)
    no = exp;
  nexp = stage ? exp : sum;
  class = 1;
  for (i = 1; i <= sum; i++)
    if (wt[i] > class)
      class = wt[i];
  fprintf(op, "%3d %5d %4d %5d %3d%3d\n", prime, nexp, facexp, no, class, 0);
  for (i = 1; i <= sum; i++)
    fprintf(op, " %3d", wt[i]);
  fprintf(op, "\n");
  if (stage == 0) {
    for (i = 1; i <= sum; i++)
      fprintf(op, " %3d", dpth[i]);
    fprintf(op, "\n");
  }
  for (i = 1; i <= sum; i++)
    fprintf(op, " %3d", d1[i]);
  fprintf(op, "\n");
  for (i = 1; i <= sum; i++)
    fprintf(op, " %3d", d2[i]);
  fprintf(op, "\n");
  jtl = 0;
  pcp = pcptr;
  e = NULL;
  while (pcp <= pcb) {
    b = *pcp;
    c = tails ? *(pcp + 1) : 0;
    if (b != 0) {
      k = *(*pcp - 1);
      l = **pcp;
      e = b + l;
    }
    else {
      k = 0;
      l = 0;
    }
    if (c != 0) {
      l += *c;
      f = c + *c;
    }
    if (k == 0 && l == 0)
      fprintf(op, "%d %d", 0, 0);
    else
      fprintf(op, "%3d  %3d", k, l);
    if (b != 0)
      while (++b <= e)
        fprintf(op, " %3d", *b);
    if (c != 0)
      while (++c <= f) {
        fprintf(op, " %3d", *c + exp);
        fprintf(op, " %3d", *(++c));
      }
    if (jtl) {
      for (i = 1; i <= 10 - l; i++)
        fprintf(op, "    ");
      if (l == 0)
        fprintf(op, "     ");
      pcpj++;
      b = *pcpj;
      if (b != 0) {
        l = *b;
        k = *(b - 1);
      }
      else {
        k = 0;
        l = 0;
      }
      e = b + l;
      fprintf(op, "%3d  %3d", k, l);
      while (++b <= e)
        fprintf(op, " %3d", *b);
    }
    fprintf(op, "\n");
    if (stage && pcp == opcb) {
      fprintf(op, " %3d\n", onng);
      for (i = 1; i <= onng; i++)
        fprintf(op, " %3d", nd1[i]);
      fprintf(op, "\n");
      for (i = 1; i <= onng; i++)
        fprintf(op, " %3d", nd2[i]);
      fprintf(op, "\n");
      if (stage == 4) {
        jump = dim * exp;
        pcpj = pcp + jump;
        jtl = 1;
      }
    }
    pcp += (1 + tails);
  }
  if (stage == 2) {
    fprintf(op, "%3d\n", chpdim);
    for (i = 1; i <= onng; i++)
      if ((b = extno[i]) != 0) {
        fprintf(op, "%3d   %4d", i, *b);
        c = b + *b;
        while (++b <= c)
          fprintf(op, " %3d", *b);
        fprintf(op, "\n");
      }
    fprintf(op, "%3d\n", chsdim);
    if (chsdim != 0)
      for (i = 1; i <= onng; i++)
        if ((b = subno[i]) != 0) {
          fprintf(op, "%3d   %4d", i, *b);
          c = b + *b;
          while (++b <= c)
            fprintf(op, " %3d", *b);
          fprintf(op, "\n");
        }
  }
  fclose(op);
  return (0);
}

void zero(short * p1, short * p2)
{
  while ((++p1) <= p2)
    *p1 = 0;
}

void setnr(short * p)
/* Really a subprocedure of collect. Adjusts nexpnt by fac times string p. */
{
  short *p1, *p2;
  p1 = p + *p;
  while (++p < p1) {
    p2 = nexpnt + *p;
    *p2 += (fac * *(++p));
    *p2 %= prime;
  }
}

int collect(short * spc, short * spf, int sgn)
/* The basic collection routine. Taken basically from Canberra NQA.
   There are several complications.
   1) The basic word is collected into an exponent vector expnt[1-exp],
      but collection also generates terms in the new generators, 1-nng.
      These are all central, and these terms come from the tails in the PCP.
      They are collected in the vector nexpnt. The procedure setnr is
      always used for this.
   2) Negative powers can occur in the string to be collected.
   3) Generators i with i>exp are generators of d(M). These do not need to be
      collected themselves; they are only there to generatoe terms in nng.
      They are not accumulated in expnt, but dealt with in blocks, as they
      occur.
   4) Generators i of P with i>facexp form an elementary abelian group. It
      saves time (particularly for larger primes) to use this fact when
      collecting
*/
{
  short gen, ex, pgen, pex, pugen, stkp, exsgn, *e, *p1, *p2, *p3, *p4, **dp,
      i, j, fpg, fpex, *p0, *pe, *q;
  char done;
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
  if (gen > exp) {
    p4 = prvec + exp + dim;
    pe = prvec + exp;
    zero(pe, p4);
    done = 0;
    exsgn = (spc < spf) ? -1 : 1;
    /* First we collect all gens i>exp in this block as a vector. We use prvec
       for this purpose, since it is not otherwise in use at present.
       (prvec is fundamentally there for use in prnrel() .)
    */
    while (1) {
      prvec[gen] += ex;
      if (spc == spf) {
        done = 1;
        break;
      }
      spc -= (2 * exsgn);
      gen = *spc;
      ex = exsgn * *(spc + 1);
      if (gen <= exp)
        break;
    }
    fpg = pgen;
    fpex = pex;
    /* Now we start to generate commutator terms. Those gens i<=facexp may
       generate new terms in prvec, so we must proceed in the correct order.
    */
    while (fpg <= facexp) {
      if (fpex > 0)
        for (i = 1; i <= fpex; i++) {
          p3 = p4;
          while (p3 > pe) {
            if ((j = *p3) != 0) {
              dp = comptr[p3 - prvec] + fpg;
              p1 = *dp;
              dp = powptr[p3 - prvec] + fpg;
              p2 = *dp;
              if (p2 != 0 && *p2 != 0) {
                fac = sgn * j;
                setnr(p2);
              }
              if (p1 != 0 && *p1 != 0) {
                p0 = p1 + 1;
                q = p1 + *p1;
                while (p0 < q) {
                  prvec[*p0] += (j * *(p0 + 1));
                  prvec[*p0] %= prime;
                  p0 += 2;
                }
              }
            }
            p3--;
          }
        }
      fpg++;
      fpex = expnt[fpg];
    }
    while (fpg <= exp) {
      if (fpex > 0) {
        p3 = p4;
        while (p3 > pe) {
          if ((j = *p3) != 0) {
            dp = comptr[p3 - prvec] + fpg;
            p1 = *dp;
            if (p1 != 0 && *p1 != 0) {
              fac = sgn * fpex * j;
              setnr(p1);
            }
          }
          p3--;
        }
      }
      fpg++;
      fpex = expnt[fpg];
    }
    p1 = expnt + exp;
    while (++pe <= p4) {
      p1++;
      *p1 += *pe;
      *p1 %= prime;
      if (*p1 < 0)
        *p1 += prime;
    }
    if (done) {
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
    sgen[stkp] = gen;
    sex[stkp] = ex;
    sspc[stkp] = spc;
  } /* gen>exp */
  else if (gen > facexp) {
    p4 = prvec + exp;
    pe = prvec + facexp;
    zero(pe, p4);
    done = 0;
    exsgn = (spc < spf) ? -1 : 1;
    /* For facexp<gen<=exp, we again treat generators in blocks, and first add
       them all up together in prvec.
    */
    while (1) {
      prvec[gen] += ex;
      if (spc == spf) {
        done = 1;
        break;
      }
      spc -= (2 * exsgn);
      gen = *spc;
      ex = exsgn * *(spc + 1);
      if (gen <= facexp || gen > exp)
        break;
    }
    fpg = pgen;
    fpex = pex;
    while (fpg <= facexp) {
      if (fpex > 0)
        for (i = 1; i <= fpex; i++) {
          p3 = tails ? p4 : prvec + no;
          while (p3 > pe) {
            if ((j = *p3) != 0) {
              dp = comptr[p3 - prvec] + (1 + tails) * fpg;
              p1 = *dp;
              p2 = *(dp + 1);
              if (tails && p2 != 0) {
                fac = sgn * j;
                setnr(p2);
              }
              if (p1 != 0 && *p1 != 0) {
                p0 = p1 + 1;
                q = p1 + *p1;
                while (p0 < q) {
                  prvec[*p0] += (j * *(p0 + 1));
                  prvec[*p0] %= prime;
                  p0 += 2;
                }
              }
            }
            p3--;
          }
        }
      fpg++;
      fpex = expnt[fpg];
    }
    /* Having generated the commutators from this block, we add it on to
       the exponent vector expnt.
    */
    p1 = expnt + facexp;
    while (++pe <= p4) {
      p1++;
      *p1 += *pe;
      *p1 %= prime;
      if (*p1 < 0)
        *p1 += prime;
    }
    if (done) {
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
    sgen[stkp] = gen;
    sex[stkp] = ex;
    sspc[stkp] = spc;
  } /* gen>facexp */

  /* Most of the remainder of the procedure follows the collect in the NQA */
  e = expnt + pgen;
  if (pgen == gen) {
    if (pugen == gen) {
      *e += ex;
      ex = 1;
      sex[stkp] = 1;
    }
    else
      (*e)++;
    if (*e < 0)
      if (gen > facexp)
        *e += prime;
      else
      /* Deal with negative exponent */
      {
        *e += prime;
        dp = powptr[gen];
        p1 = *dp;
        p2 = *(dp + 1);
        if (tails && p2 != 0) {
          fac = -sgn;
          setnr(p2);
        }
        if (p1 != 0 && *p1 != 0) {
          exsgn = -1;
          pex = -1;
          spc = p1 + 1;
          spf = spc + *p1 - 2;
          goto recurse;
        }
      }
    else if (*e >= prime) {
      if (gen > facexp)
        *e -= prime;
      else
      /* Deal with exponent >= prime */
      {
        *e -= prime;
        dp = powptr[gen];
        p1 = *dp;
        p2 = *(dp + 1);
        if (tails && p2 != 0) {
          fac = sgn;
          setnr(p2);
        }
        if (p1 != 0 && *p1 != 0) {
          exsgn = 1;
          pex = -1;
          pugen = gen;
          spf = p1 + 1;
          spc = spf + *p1 - 2;
          goto recurse;
        }
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
  else if (pgen > facexp) {
    pex = 0;
    goto loop;
  }
  else {
    dp = comptr[gen] + (1 + tails) * pgen;
    p1 = *dp;
    p2 = *(dp + 1);
    if (p1 != 0 && *p1 != 0) {
      if (pugen == pgen) {
        spugen[stkp] = -pgen;
        if (ex < 0) {
          if (gen > facexp)
            sex[stkp] = ex + prime;
          else {
            sex[stkp] = ex + prime;
            dp = powptr[gen];
            p3 = *dp;
            p4 = *(dp + 1);
            if (tails && p4 != 0) {
              fac = -sgn;
              setnr(p4);
            }
            if (p3 != 0 && *p3 != 0) {
              exsgn = -1;
              spc = p3 + 1;
              spf = spc + *p3 - 2;
              goto recurse;
            }
          }
        }
      }
      else
        pugen = pgen;
    }
    else {
      if (tails && p2 != 0) {
        fac = (pugen == pgen) ? ex * pex * sgn : pex * sgn;
        setnr(p2);
      }
      pex = 0;
      goto loop;
    }
    spf = p1 + 1;
    spc = spf + *p1 - 2;
    exsgn = 1;
    if (tails && p2 != 0) {
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

void setpinv(void)
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
}

void bgc(void)
/* Garbage collection by outputting to a temporary file and reading in again
   Occurs after wasted space accumulates from eliminating generators.
*/
{
  short **pcp, *p, *q, **epcp;
  int     l, k, d, ct;
  char    gcfile[80];
  FILE *  wt;
  strcpy(gcfile, inf);
  strcat(gcfile, "wrd.t");
  printf("Back garbage collection. wsp=%d.\n", wsp);
  fflush(stdout);
  fflush(stderr);
  wt = fopen(gcfile, "w");
  if (stage >= 2) {
    pcp = pcb + 1;
    ct = 0;
    epcp = (stage > 2) ? npcb2 : npcb;
    d = 1;
  }
  else if (stage == 1) {
    pcp = opcb + facexp + 1;
    ct = facexp + 1;
    epcp = pcb;
    d = 1;
  }
  else {
    pcp = pcptr;
    ct = 0;
    epcp = pcb;
    d = 0;
  }
  while (pcp <= epcp) {
    if (stage == 0)
      pcp++;
    if ((p = *pcp) != 0) {
      q = p + *p;
      p -= (1 + d);
      while (++p <= q)
        fprintf(wt, "%d ", *p);
    }
    pcp++;
    if (stage == 1) {
      if (ct == exp) {
        pcp += facexp;
        ct = facexp + 1;
      }
      else
        ct++;
    }
  }
  rpb = rel + rsp - 1;
  printf("Written to file.\n");
  fflush(stdout);
  fclose(wt);
  wt = fopen(gcfile, "r");
  if (stage >= 2) {
    pcp = pcb + 1;
    epcp = (stage > 2) ? npcb2 : npcb;
  }
  else if (stage == 1) {
    pcp = opcb + facexp + 1;
    ct = facexp + 1;
    epcp = pcb;
  }
  else {
    pcp = pcptr;
    epcp = pcb;
  }
  while (pcp <= epcp) {
    if (stage == 0)
      pcp++;
    if (*pcp != 0) {
      if (stage) {
        fscanf(wt, "%d%d", &k, &l);
        rpb -= (l + 2);
        p = rpb + 2;
        *(p - 1) = k;
      }
      else {
        fscanf(wt, "%d", &l);
        rpb -= (l + 1);
        p = rpb + 1;
      }
      *pcp = p;
      *p = l;
      q = p + l;
      while (++p <= q)
        fscanf(wt, "%hd", p);
    }
    pcp++;
    if (stage == 1) {
      if (ct == exp) {
        pcp += facexp;
        ct = facexp + 1;
      }
      else
        ct++;
    }
  }
  wsp = 0;
  fclose(wt);
  unlink(gcfile);
  printf("Reread file.\n");
  fflush(stdout);
}

int intgen(int i, int j)
/* A new generator is introduced for commutator [i,j] or (if i=j) i^p */
{
  short ** dp;
  int      sum;
  dp = (stage) ? comptr[i] + j : (i == j) ? powptr[i] : comptr[i] + 2 * j;
  if (*dp != 0 &&
      *(*dp - 1) != 0) { /*printf("Reln [%d,%d] is already a defn.\n",i,j);*/
    return (1);
  }
  if (stage >= 2)
    dp = powptr[i] + j;
  if (stage == 0)
    dp++;
  if (*dp !=
      0) { /*fprintf(stderr,"Reln [%d,%d] already has a tail.\n",i,j);*/
    return (1);
  }
  nng++;
  enexpnt++;
  if (nng >= mng) {
    fprintf(stderr, "Too many new gens. Increase MNG.\n");
    return (-1);
  }
  if (stage) {
    nd1[nng] = i;
    nd2[nng] = j;
    rpb -= 4;
    *(rpb + 1) = nng;
    *(rpb + 2) = 2;
    *(rpb + 3) = nng;
    *(rpb + 4) = 1;
    *dp = rpb + 2;
  }
  else {
    sum = exp + nng;
    d1[sum] = i;
    d2[sum] = j;
    dpth[nng + exp] = depth;
    wt[sum] = (i == j) ? wt[i] + 1 : wt[i] + wt[j];
    rpb -= 3;
    *(rpb + 1) = 2;
    *(rpb + 2) = nng;
    *(rpb + 3) = 1;
    *dp = rpb + 1;
    dp--;
    if (*dp == 0) {
      rpf++;
      *dp = rpf;
      *rpf = 0;
      rpf++;
    }
    *(*dp - 1) = sum;
  }
  if (rpb + 1 - rpf < marg) {
    printf("Running out of space in intgen. Remains,marg,wsp=%d,%d,%d.\n",
           (int)(rpb + 1 - rpf), marg, wsp);
    if (wsp > marg)
      bgc();
    else
      return (-1);
  }
  return (0);
}

int subrel(int i, int j)
/* Value of [i,j] or i^p is computed using an associativity condition. */
{
  short ** dp;
  int      x, y, z;
  dp = (stage) ? comptr[i] + j : (i == j) ? powptr[i] : comptr[i] + 2 * j;
  if (*dp != 0 &&
      *(*dp - 1) != 0) { /*printf("Reln [%d,%d] is already a defn.\n",i,j);*/
    return (0);
  }
  if (stage >= 2)
    dp = powptr[i] + j;
  if (stage == 0)
    dp++;
  if (*dp != 0) {
    fprintf(stderr, "Reln [%d,%d] already has a tail.\n", i, j);
    return (-1);
  }
  z = 0;
  for (x = nng; x >= 1; x--)
    if ((y = nexpnt[x]) != 0) {
      if (y < 0)
        y += prime;
      *rpb = y;
      *(rpb - 1) = x;
      rpb -= 2;
      z += 2;
    }
  *rpb = z;
  if (stage) {
    *dp = rpb;
    rpb -= 2;
    *(rpb + 1) = 0;
  }
  else {
    *dp = rpb;
    rpb--;
  }
  if (rpb + 1 - rpf < marg) {
    printf("Running out of space in subrel. Remains,marg,wsp=%d,%d,%d.\n",
           (int)(rpb + 1 - rpf), marg, wsp);
    if (wsp > marg)
      bgc();
    else
      return (-1);
  }
  return (0);
}

int assoc(int g1, int g2, int g3)
/* Associativity relation is collected. */
{
  char    eq12, eq23, prnt, triv;
  short * p;
  int     i, e, sum;
  sum = stage ? exp + dim : exp;
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
  for (i = 1; i <= sum; i++) {
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

int prnrel(void)
/* The new relation in nexpnt is processed, and a new generator eliminated
   or made redundant if necessary.
*/
{
  int    i, j, k, l, nl, w, x, y, ct, elno, fac, pow, len, gno;
  short *p, *q, *eprvec, **dp, **edp;
  char   sub, gth, u;
  u = (stage) ? 1 : 0;
  pow = 0;
  len = 0;
  gno = 0;
  for (i = nng; i >= 1; i--)
    if ((x = nexpnt[i]) != 0) {
      len++;
      if (pow == 0 || (stage == 0 && wt[i + exp] < w)) {
        pow = x;
        gno = i;
        w = wt[i + exp];
      }
    }
  if (len == 0)
    return (1);
  /* When stage=2 and gno<=onng, we do not eliminate a generator. In this
     case, the relation is in fact a generator of the dual of H^2(P,M), and it
     is stored as such. Otherwise we will be eliminating gno.
  */
  if (stage == 2 && gno <= onng)
    return (1 + gno);
  sub = len > 1;
  if (sub) {
    if (pow < 0)
      pow += prime;
    fac = pinv[pow];
    p = nexpnt;
    while (++p <= enexpnt) {
      *p *= fac;
      *p %= prime;
    }
    /* printf("New gen no %d is redundant.\n",gno); */
  }
  /* else printf("New gen no %d is eliminated.\n",gno); */
  eprvec = prvec + nng;
  elno = (u) ? gno : exp + gno;
  if (stage >= 2) {
    dp = pcb + 1;
    ct = 0;
    edp = (stage > 2) ? npcb2 : npcb;
  }
  else if (stage == 1) {
    dp = opcb + facexp + 1;
    ct = facexp + 1;
    edp = pcb;
  }
  else {
    dp = pcptr;
    ct = 0;
    edp = pcb;
  }
  /* Now we edit the PCP relations, wasting as little space as possible */
  while (dp <= edp) {
    p = *dp;
    if (p != 0) {
      p--;
      if (*p > elno)
        (*p)--;
      else if (*p == elno)
        *p = 0;
      p++;
    }
    if (stage == 0) {
      dp++;
      p = *dp;
    }
    if (p != 0) {
      l = *p;
      q = p + l;
      gth = 0;
      while (++p < q) {
        if (gth) {
          *(p - 2) = *p - 1;
          *(p - 1) = *(p + 1);
        }
        else if (*p > gno)
          (*p)--;
        else if (*p == gno) {
          gth = 1;
          if (sub)
            break;
          (**dp) -= 2;
          if (**dp == 0) {
            *dp = 0;
            wsp += (3 + u);
          }
          else
            wsp += (2 + u);
        }
        p++;
      }
      if (gth && sub) {
        zero(prvec, eprvec);
        p = *dp;
        while (++p < q) {
          prvec[*p] = *(p + 1);
          p++;
        }
        y = prvec[gno];
        nl = l - 2;
        for (j = 1; j <= nng; j++)
          if (j != gno && (k = nexpnt[j]) != 0) {
            q = prvec + j;
            if (*q == 0)
              nl += 2;
            *q -= y * k;
            *q %= prime;
            if (*q < 0)
              *q += prime;
            if (*q == 0)
              nl -= 2;
          }
        if (nl == 0) {
          *dp = 0;
          wsp += (l + 1 + u);
        }
        else {
          if (nl > l) {
            if (rpb + 1 - rpf < marg) {
              printf("Running out of space in prnrel. "
                     "Remains,marg,wsp=%d,%d,%d.\n",
                     (int)(rpb + 1 - rpf), marg, wsp);
              if (wsp > marg)
                bgc();
              else
                return (-1);
            }
            if (stage) {
              rpb -= (nl + 2);
              *(rpb + 1) = *(*dp - 1);
              *dp = rpb + 2;
              wsp += (l + 2);
            }
            else {
              rpb -= (nl + 1);
              *dp = rpb + 1;
              wsp += (l + 1);
            }
          }
          for (j = gno; j < nng; j++)
            prvec[j] = prvec[j + 1];
          p = *dp;
          *p = nl;
          for (j = 1; j < nng; j++)
            if ((y = prvec[j]) > 0) {
              *(++p) = j;
              *(++p) = y;
            }
        }
      }
    }
    dp++;
    if (stage == 1) {
      if (ct == exp) {
        dp += facexp;
        ct = facexp + 1;
      }
      else
        ct++;
    }
  }
  for (i = gno; i < nng; i++)
    if (stage) {
      nd1[i] = nd1[i + 1];
      nd2[i] = nd2[i + 1];
    }
    else {
      j = exp + i;
      k = j + 1;
      d1[j] = d1[k];
      d2[j] = d2[k];
      wt[j] = wt[k];
    }
  nng--;
  enexpnt--;
  if (nng == 0) {
    printf("All new generators eliminated.\n");
    return (0);
  }
  return (1);
}
