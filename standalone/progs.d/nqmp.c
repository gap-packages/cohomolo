#include "defs.h"

extern char  ims, act, gap, crel, inf0[], inf1[], inf2[], outf[], outfm[];
extern short intexp, mng, mexp, wksp, prime, exp, nng, class, *rpf, *rpb,
    *eexpnt, *enexpnt, **pcb, mnng, mord, rel[], expnt[], nexpnt[], cord[],
    wt[], d1[], d2[], *pcptr[], **powptr[], **comptr[], *sspc[], *sspf[],
    sgen[], sex[], spgen[], spex[], spugen[], *intg[], *imintg[], *tlintg[];
extern int ptrsp, rsp;
short *    wf, *wc;
char       norm;

/* The data structures for this program and for nqrun are similar.
   d1 and d2 contain definitions of generators. (Def. comes from commutator
   [d1[i],d2[i]], or if d1[i]=d2[i] from the p-th power of d1[i].)
   The PCP is stored in rel. Commutators and powers are pointed to by comptr
   and powptr. More precisely, they are pointed to by pointers stored in
   pcptr, and it is these pointers that are pointed to by comptr and powptr.
   Usually, each such relation has a tail, which is a word in the
   new generators. If so, comptr[i]+2*j and comptr[i]+2*j+1 point to the
   relation in P, and its tail, respectively, and similarly for powptr.
   Relations are stored as gen-pow strings, preceded by length. (The length
   entry is what is pointed to.) This is also preceded by a number
   indicating if this is a definition. In nqmrun, the tails are stored as
   vectors in the new generators, but in nqrun, they are stored as gen-pow
   strings, to save space.
   When collecting, the collected element is computed into the array expnt, as
   a vector, and its tail is computed into nexpnt. eexpnt and enexpnt point to
   the ends of these vectors.
   Details of file formats are in Info.5
*/
int nqmprog(void)
{
  short i, j, k, l, m, *gi, *gj, *ti, *tj, cl, def, *ps, *pf, **dp, *nrpb, *p,
      *orpf, *orpb, nb, np, k1, *rno, *covrel, **pgen, sgn;
  char  nt;
  FILE *ip, *op;
  if (ingp() == -1) {
    fprintf(stderr, "Input error.\n");
    return (-1);
  }
  covrel = NULL;
  eexpnt = expnt + exp;
  enexpnt = nexpnt + nng;

  /* if nng=0, we are computing the multiplier of P. First we estimate the
     maximal possible no of new gens, mnng.
  */
  if (nng == 0) {
    mnng = 0;
    for (i = 1; i <= exp; i++)
      if (wt[i] > 1)
        mnng++;
    for (i = 1; i <= exp; i++)
      for (j = 1; j < i; j++) {
        gi = *(comptr[i] + 2 * j);
        if ((wt[i] + wt[j] <= class + 1) && (wt[i] == 1 || wt[j] == 1) &&
            (gi == 0 || *(gi - 1) == 0))
          mnng++;
      }
    printf("mnng=%d\n", mnng);
    j = exp / 2;
    mord = 1;
    for (i = 1; i <= j; i++)
      mord *= prime;
    /* Now we start introducing the new generators */
    for (cl = class + 1; cl >= 2; cl--) {
      for (i = 1; i < exp; i++)
        for (j = i + 1; j <= exp; j++)
          if (wt[i] + wt[j] == cl && (wt[i] == 1 || wt[j] == 1))
            if (intgen(j, i) == -1)
              return (-1);
      for (i = 3; i <= exp; i++)
        if (wt[i] == cl)
          if (intgen(i, i) == -1)
            return (-1);
      for (i = 3; i < exp; i++)
        if (wt[i] > 1)
          for (j = i + 1; j <= exp; j++)
            if (wt[j] > 1 && wt[i] + wt[j] == cl) {
              k = d1[i];
              l = d2[i];
              if (assoc(j, k, l))
                if (subrel(j, i) == -1)
                  return (-1);
            }
      for (j = 1; j <= exp; j++)
        for (i = 1; i <= j; i++)
          if (wt[i] + wt[j] == cl) {
            if (assoc(j, i, i))
              if ((l = prnrel(0)) == -1)
                goto nextcl;
            if (j != i && assoc(j, j, i))
              if ((l = prnrel(0)) == -1)
                goto nextcl;
          }
      for (k = 3; k <= exp; k++)
        for (j = 2; j < k; j++)
          for (i = 1; i < j; i++)
            if (wt[i] + wt[j] + wt[k] == cl && assoc(k, j, i))
              if ((l = prnrel(0)) == -1)
                goto nextcl;
    nextcl:;
    }
    if (nng == 0) {
      if (gap) {
        op = fopen(outfm, "w");
        fprintf(op, "COHOMOLO.Multiplier:=[];\n");
        fclose(op);
        printf("Multiplier is trivial.\n");
      }
      else
        fprintf(stderr, "Multiplier is trivial.\n");
      return (-1);
    }
    /* Multiplier of P is now computed */
  }
  if (act) {
    orpf = rpf;
    orpb = rpb;
    if ((ip = fopen(inf2, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf2);
      return (-1);
    }
    fscanf(ip, "%hd", &intexp);
    while (intexp != -1) {
      norm = intexp == exp;
      if (ims) {
        if (norm == 0) {
          fprintf(stderr, "-i can only be called for g in N(P).\n");
          return (-1);
        }
        printf("\nAction of g on multiplier is as follows:\n\n");
      }
      /* First we read in the genrators and their images under the action.
         These images will later have their own tail in the new generators,
         which will be stored in the arrays tlintg[i].
      */
      rpf = orpf;
      rpb = orpb;
      if (norm == 0)
        for (i = 1; i <= intexp; i++) {
          fscanf(ip, "%hd", rpf);
          intg[i] = rpf;
          l = *rpf;
          for (j = 1; j <= l; j++) {
            rpf++;
            fscanf(ip, "%hd", rpf);
          }
          rpf++;
        }
      for (i = 1; i <= intexp; i++) {
        fscanf(ip, "%hd", rpf);
        imintg[i] = rpf;
        tlintg[i] = 0;
        l = *rpf;
        for (j = 1; j <= l; j++) {
          rpf++;
          fscanf(ip, "%hd", rpf);
        }
        rpf++;
      }
      if (rpb - rpf < wksp) {
        fprintf(stderr, "Out of space. Increase RSP.\n");
        return (-1);
      }
      /* Now we compute the action */
      if (norm) {
        wf = rpf;
        for (l = 3; l <= exp; l++)
          if ((i = d1[l]) != 0) {
            j = d2[l];
            p = *(comptr[i] + 2 * j);
            ps = p + 1;
            pf = p + *p - 1;
            zero(nexpnt, enexpnt);
            p = pf - 2;
            wc = wf - 2;
            while (p >= ps) {
              enter(imintg[*p], -*(p + 1), tlintg[*p]);
              p -= 2;
            }
            gi = imintg[i];
            gj = imintg[j];
            ti = tlintg[i];
            tj = tlintg[j];
            enter(gi, -1, ti);
            enter(gj, -1, tj);
            enter(gi, 1, ti);
            enter(gj, 1, tj);
            zero(expnt, eexpnt);
            collect(wc, wf, 1);
            nt = 0;
            for (k = 1; k <= nng; k++)
              if (nexpnt[k] != 0) {
                nt = 1;
                break;
              }
            if (nt) {
              nrpb = rpb - nng;
              zero(nrpb, rpb);
              tlintg[l] = nrpb;
              for (k = 1; k <= nng; k++)
                nrpb[k] = nexpnt[k];
              rpb = nrpb;
            }
          }
      }
      for (i = 1; i <= intexp; i++)
        for (j = 1; j <= i; j++) {
          if (ims == 0 && i == j)
            continue;
          if (norm) {
            if (ims && i == j)
              dp = powptr[i];
            else
              dp = comptr[i] + 2 * j;
            p = *dp;
            if (p == 0) {
              def = 0;
              ps = rpf;
              pf = ps - 2;
            }
            else {
              if ((def = *(p - 1)))
                continue;
              ps = p + 1;
              pf = p + *p - 1;
            }
            covrel = *(dp + 1);
            if (ims)
              if (covrel == 0 || *covrel == 0)
                continue;
          }
          else {
            if (i == intexp) {
              def = 0;
              *rpf = 0;
            }
            else {
              fscanf(ip, "%hd%hd", &def, rpf);
              for (k = 1; k <= *rpf; k++)
                fscanf(ip, "%hd", rpf + k);
            }
            ps = rpf + 1;
            wf = ps + *rpf;
            pf = wf - 2;
          }
          if (def) {
            zero(nexpnt, enexpnt);
            p = pf - 2;
            wc = wf - 2;
            while (p >= ps) {
              enter(imintg[*p], -*(p + 1), tlintg[*p]);
              p -= 2;
            }
            gi = imintg[i];
            gj = imintg[j];
            ti = tlintg[i];
            tj = tlintg[j];
            enter(gi, -1, ti);
            enter(gj, -1, tj);
            enter(gi, 1, ti);
            enter(gj, 1, tj);
            zero(expnt, eexpnt);
            collect(wc, wf, 1);
            p = pf - 2;
            wc = wf - 2;
            while (p >= ps) {
              enter(intg[*p], -*(p + 1), 0);
              p -= 2;
            }
            gi = intg[i];
            gj = intg[j];
            enter(gi, -1, 0);
            enter(gj, -1, 0);
            enter(gi, 1, 0);
            enter(gj, 1, 0);
            zero(expnt, eexpnt);
            collect(wc, wf, -1);
            nt = 0;
            for (k = 1; k <= nng; k++)
              if (nexpnt[k] != 0) {
                nt = 1;
                break;
              }
            if (nt) {
              nrpb = rpb - nng;
              zero(nrpb, rpb);
              tlintg[def] = nrpb;
              for (k = 1; k <= nng; k++)
                nrpb[k] = nexpnt[k];
              rpb = nrpb;
            }
          }
          else {
            zero(nexpnt, enexpnt);
            zero(expnt, eexpnt);
            if (norm == 0) {
              p = pf;
              wc = wf - 2;
              while (p >= ps) {
                enter(intg[*p], -*(p + 1), 0);
                p -= 2;
              }
              gi = intg[i];
              gj = intg[j];
              enter(gi, -1, 0);
              enter(gj, -1, 0);
              enter(gi, 1, 0);
              enter(gj, 1, 0);
              collect(wc, wf, -1);
            }
            p = pf;
            wc = wf - 2;
            while (p >= ps) {
              enter(imintg[*p], -*(p + 1), tlintg[*p]);
              p -= 2;
            }
            if (i == j)
              enter(imintg[i], prime, tlintg[i]);
            else {
              gi = imintg[i];
              gj = imintg[j];
              ti = tlintg[i];
              tj = tlintg[j];
              enter(gi, -1, ti);
              enter(gj, -1, tj);
              enter(gi, 1, ti);
              enter(gj, 1, tj);
            }
            zero(expnt, eexpnt);
            collect(wc, wf, 1);
            if (ims) {
              for (k = 1; k <= nng; k++)
                printf("%2d", covrel[k]);
              printf("    ->    ");
              for (k = 1; k <= nng; k++) {
                m = nexpnt[k];
                if (m < 0) {
                  l = cord[k];
                  m += l;
                }
                printf("%2d", m);
              }
              printf("\n");
            }
            else {
              if (norm && covrel)
                for (k = 1; k <= nng; k++) {
                  l = covrel[k];
                  nexpnt[k] -= l;
                  l = cord[k];
                  nexpnt[k] %= l;
                }
              if ((l = prnrel(1)) == -1)
                return (-1);
            }
          }
        }
      fscanf(ip, "%hd", &intexp);
    }
    fclose(ip);
  }
  outgp();
  if (gap) {
    op = fopen(outfm, "w");
    fprintf(op, "COHOMOLO.Multiplier:=[");
    for (i = 1; i <= nng; i++) {
      fprintf(op, "%d", cord[i]);
      if (i < nng)
        putc(',', op);
    }
    fprintf(op, "];\n");
    fclose(op);
    printf("Orders of the cyclic factors of the multiplier are now:\n");
    for (i = 1; i <= nng; i++)
      printf("%3d  ", cord[i]);
    printf("\n");
  }
  else {
    fprintf(stderr,
            "Orders of the cyclic factors of the multiplier are now:\n");
    for (i = 1; i <= nng; i++)
      fprintf(stderr, "%3d  ", cord[i]);
    fprintf(stderr, "\n");
  }
  if (crel) {
    strcpy(inf1, inf0);
    strcat(inf1, "psgwds");
    if ((ip = fopen(inf1, "r")) == 0) {
      printf("Cannot open %s.\n", inf1);
      return (-1);
    }
    pgen = pcb;
    fscanf(ip, "%hd", &np);
    if (pcb + np - pcptr >= ptrsp) {
      fprintf(stderr, "Out of ptr space. Increase PTRSP.\n");
      return (-1);
    }
    for (i = 1; i <= np; i++) {
      pgen[i] = rpf;
      fscanf(ip, "%hd", rpf);
      p = rpf;
      rpf += (1 + *p);
      while (++p < rpf)
        fscanf(ip, "%hd", p);
    }
    fclose(ip);
    strcpy(inf1, inf0);
    strcat(inf1, "psg.rel");
    if ((ip = fopen(inf1, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf1);
      return (-1);
    }
    strcpy(outf, inf0);
    strcat(outf, "psg.er");
    op = fopen(outf, "w");
    fscanf(ip, "%hd", &nb);
    rno = rpf - 1;
    rpf += nb;
    if (rpb - rpf < wksp) {
      fprintf(stderr, "Out of space. Increase RSP.\n");
      return (-1);
    }
    for (i = 1; i <= nb; i++)
      fscanf(ip, "%hd", rno + i);
    fprintf(op, "%4d%4d\n", nb, nng);
    for (i = 1; i <= nb; i++)
      fprintf(op, "%4d", rno[i]);
    fprintf(op, "\n");
    for (i = 1; i <= nng; i++)
      fprintf(op, "%4d", cord[i]);
    fprintf(op, "\n");
    wf = rpf;
    for (i = 1; i <= rno[1]; i++) {
      fscanf(ip, "%hd", &l);
      covrel = rpb - l;
      p = covrel;
      while (++p <= rpb)
        fscanf(ip, "%hd", p);
      zero(expnt, eexpnt);
      zero(nexpnt, enexpnt);
      for (j = l; j >= 1; j--) {
        wc = wf - 2;
        k = covrel[j];
        k1 = k / 2 + 1;
        m = *pgen[k1];
        if (k % 2 == 0) {
          sgn = 1;
          ps = pgen[k1] + 1;
          pf = ps + m - 2;
        }
        else {
          sgn = -1;
          pf = pgen[k1] + 1;
          ps = pf + m - 2;
        }
        while (1) {
          wc += 2;
          *wc = *ps;
          *(wc + 1) = *(ps + 1) * sgn;
          if (ps == pf)
            break;
          ps += (2 * sgn);
        }
        collect(wc, wf, 1);
      }
      fprintf(op, "%4d  ", l);
      for (j = 1; j <= l; j++)
        fprintf(op, "%4d", covrel[j]);
      fprintf(op, "\n");
      l = 0;
      p = nexpnt;
      while (++p <= enexpnt)
        if (*p != 0) {
          l += 2;
          if (*p < 0)
            *p += cord[p - nexpnt];
        }
      fprintf(op, "%4d  ", l);
      p = nexpnt;
      while (++p <= enexpnt)
        if (*p != 0)
          fprintf(op, "%4d%4d", (int)(p - nexpnt), *p);
      fprintf(op, "\n");
      if ((i == rno[1]) && (fscanf(ip, "%hd", &j) > 0)) {
        fprintf(op, "%4d\n", j);
        rno[1] += j;
      }
    }
  }
  return (0);
}

int enter(short * g, int pow, short * t)
/* This enters a power of the gen-pow string pointed to by g into the word
   which we are building up for collection.
*/
{
  short *ps, *pf, *pc, i, j, sgn;
  if (t != 0)
    for (i = 1; i <= nng; i++) {
      nexpnt[i] += (pow * t[i]);
      j = cord[i];
      nexpnt[i] %= j;
    }
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
  return (0);
}
