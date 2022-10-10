#include "defs.h"

extern char inf1[], inf2[], inf3[], outf1[], outf2[], outf3[], temp1[],
    temp2[], expg, exph, cr, dcr, triv;
extern int psp, trsp, svsp;
extern int mpt, mb, tree[], perm[], gorb[], lorbg[], lorbh[], base[], fpt[],
    bpt[], coh_index[], *cp[], *trad[], *sv[], *tailad[], **cpad[],
    **svgptr[], **svhptr[];
int trptr;
int npt, **lcp, **ucp, *stop, dummy, *pst, *pend, nb, fnt, lnt, ind, coind,
    oind, cind, bno, pt, **expp, npg;
FILE *ip, *op, *opx;

void advance(void);

/* The data structures of this program differ from other permutation group
   programs. Permutations are not numbered, but are located by their base
   adresses. Consequently, arrays sv and cp are arrays of pointers. The
   current permutation in cp is stored in the middle of cp (since it has to be
   expanded in both directions), between adresses lcp and ucp of cp.
   The next six procedures are the equivlents of the corresponding ones
   in permfns.c, but addsv can go in both directions.
   In this half of the program, coset reps of H in G are computed.
*/

int image(int pt)
{
  int ** p;
  p = lcp;
  while (p <= ucp) {
    pt = (*p)[pt];
    p++;
  }
  return (pt);
}

void addsvb(int pt, int ** sv)
{
  int * p;
  while ((p = sv[pt]) != stop) {
    pt = p[pt];
    lcp--;
    *lcp = p - npt;
  }
}

void addsvf(int pt, int ** sv)
{
  int * p;
  while ((p = sv[pt]) != stop) {
    pt = p[pt];
    ucp++;
    *ucp = p;
  }
}

void invert(int * a, int * b)
{
  int i;
  for (i = 1; i <= npt; i++)
    b[a[i]] = i;
}

void rdperm(int * a, int * b)
{
  int i, *r;
  for (i = 1; i <= npt; i++)
    fscanf(ip, "%d", a + i);
  r = pst + 1;
  fscanf(ip, "%d", r);
  invert(a, b);
}

void rdsv(int ** sv)
{
  int i, j;
  for (i = 1; i <= npt; i++) {
    fscanf(ip, "%d", &j);
    sv[i] = (j == 0) ? 0 : (j == -1) ? stop : cp[j];
  }
}

int cnprg1(void)
{
  int  i, j, k, ad, sad, nph, *svpt, *p1, *p2, **csvg, **csvh, nexp;
  char pthere, w, ok, allok;
  long fac;
  stop = &dummy;
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%d%d%d%d", &npt, &npg, &nb, &k);
  if (npt > mpt) {
    fprintf(stderr, "npt too big. Increase mpt.\n");
    return (-1);
  }
  if (nb >= mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (2 * nb * npt > svsp) {
    fprintf(stderr, "svsp too small. Increase SVSP.\n");
    return (-1);
  }
  if (k <= 0) {
    fprintf(stderr, "Wrong input format for G.\n");
    return (-1);
  }
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%d", base + i);
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%d", lorbg + i);
  pst = perm - 1;
  pend = perm + psp - 1;
#ifdef __clang__
  /* workaround for a bug in Apple LLVM version 7.0.2 (clang-700.1.81),
    and possibly other versions. See also
    https://github.com/gap-packages/cohomolo/issues/5
  */
  if (psp < 0) /* never should happen */
    fprintf(stderr, "%d", psp);
#endif
  if (pst + 2 * npg * npt > pend) {
    fprintf(stderr, "Out of space.Increase PSP.\n");
    return (-1);
  }
  for (i = 1; i <= npg; i++) {
    p1 = pst;
    p2 = pst + npt;
    pst += (2 * npt);
    rdperm(p1, p2);
    cp[2 * i - 1] = p2;
  }
  for (i = 1; i <= nb; i++) {
    svgptr[i] = sv + npt * (i - 1) - 1;
    rdsv(svgptr[i]);
  }
  fclose(ip);
  /* G is read. Now read the subgroup H */
  if (triv) {
    nph = 0;
    for (i = 1; i <= nb; i++) {
      lorbh[i] = 1;
      svhptr[i] = sv + npt * (nb + i - 1) - 1;
      for (j = 1; j <= npt; j++)
        svhptr[i][j] = 0;
      svhptr[i][base[i]] = stop;
    }
  }
  else {
    if ((ip = fopen(inf2, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf2);
      return (-1);
    }
    fscanf(ip, "%d%d%d%d", &i, &nph, &j, &k);
    w = 0;
    if (k <= 0) {
      fprintf(stderr, "Wrong input format for H.\n");
      return (-1);
    }
    if (i != npt || j != nb)
      w = 1;
    else
      for (i = 1; i <= nb; i++) {
        fscanf(ip, "%d", lorbh + i);
        if (lorbh[i] != base[i])
          w = 1;
      }
    if (w) {
      fprintf(stderr, "npt or base for G and H do not agree.\n");
      return (-1);
    }
    for (i = 1; i <= nb; i++)
      fscanf(ip, "%d", lorbh + i);
    if (pst + 2 * npt * nph + 2 > pend) {
      fprintf(stderr, "Out of space.Increase PSP.\n");
      return (-1);
    }
    for (i = 1; i <= nph; i++) {
      p2 = pend - npt;
      p1 = p2 - npt;
      pend -= (2 * npt);
      rdperm(p1, p2);
      cp[2 * i - 1] = p2;
    }
    for (i = 1; i <= nb; i++) {
      svhptr[i] = sv + npt * (nb + i - 1) - 1;
      rdsv(svhptr[i]);
    }
    fclose(ip);
  }
  for (i = 1; i <= nb; i++)
    fpt[i] = 0;
  p1 = pst;
  p2 = p1 + npt;
  pst += (2 * npt);
  lcp = cp;
  expp = cp + 2 * npt - 1;
  cind = 1;
  fac = 1;
  coh_index[nb + 1] = 1;
  lnt = 0;

  /* Now we compute the indices of H in G in the stabilizer chain into the
     array coh_index. The nontrivial indices are linked by fpt and bpt.
  */
  for (i = nb; i >= 1; i--) {
    j = lorbg[i];
    fac = fac * j;
    j = lorbh[i];
    fac = fac / j;
    ind = fac;
    coh_index[i] = ind;
    if (ind > coh_index[i + 1]) {
      if (lnt == 0)
        lnt = i;
      else {
        fpt[i] = fnt;
        bpt[fnt] = i;
      }
      fnt = i;
    }
  }
  bpt[fnt] = -1;
  fpt[lnt] = -1;
  printf("Indices:  ");
  for (i = 1; i <= nb; i++)
    printf("%5d", coh_index[i]);
  printf("\n");
  if (ind == 1) {
    fprintf(stderr, "G=H.\n");
    return (-1);
  }
  nexp = 0;
  ind = 1;
  oind = 1;
  trptr = 0;

  /* Now we start the backtrack search for the left coset reps of H in G. The
     necessary information about these is stored in the array tree. This is
     done so as to make it as easy as possible to compute which coset an
     arbitrary element of G lies in. If expg is true, then those perms in G
     which arise are expanded and stored. If dcr is true, then a small subset
     of the coset reps will need to be output as perms, later. For this
     purpose, the base images of all coset reps are output to a temporary
     file. ind is the current coh_index coh_index[bno], and cind is the no of
     coset reps found so far.
  */
  for (i = fnt; i != -1; i = fpt[i]) {
    trad[i] = tree + trptr;
    tree[trptr] = base[i];
    trptr += 2;
  }
  if (dcr)
    opx = fopen(temp2, "w");
  if (cr) {
    op = fopen(outf3, "w");
    fprintf(op, "%4d %d%4d%4d\n", npt, coh_index[1] - 1, 0, 0);
  }
  sad = trptr;
  tree[trptr] = 1;
  trptr++;
  for (bno = lnt; bno != -1; bno = bpt[bno]) {
    printf("bno=%d.\n", bno);
    *gorb = 1;
    gorb[1] = base[bno];
    ind = coh_index[bno];
    csvg = svgptr[bno];
    csvh = svhptr[bno];
    allok = (lorbh[bno] == 1);
    /* if allok, then all coset reps in this link of the stabilizer chain will
       be coset reps of H in G.
    */
    if (expg) {
      for (i = 1; i <= npt; i++)
        expp[i] = 0;
      expp[base[bno]] = stop;
    }
    for (pt = 1; pt <= npt; pt++) {
      svpt = csvg[pt];
      if (svpt != 0 && svpt != stop) {
        ucp = lcp - 1;
        addsvf(pt, csvg);
        pthere = 0;
        if (expg) {
          for (i = 1; i <= npt; i++)
            p1[i] = image(i);
          invert(p1, p2);
          ucp = lcp;
          *ucp = p1;
        }
        j = sad;
        for (i = fpt[bno]; i != -1; i = fpt[i]) {
          tailad[i] = tree + j;
          cpad[i] = ucp + 1;
          j += 2;
        }
        /* In the search, for each basic coset rep g in the stabilizer chain,
           we have to try out gh as a possible coset rep of H in G, for each h
           in the list of coset reps of H[bno-1] in G[bno-1]. coind and oind
           record the search through this list.
        */
        coind = 1;
        while (coind <= oind) {
          if (coind > 1)
            advance();
          ok = 1;
          if (allok == 0)
            for (i = 1; i <= *gorb && ok; i++)
              if (csvh[image(gorb[i])] != 0)
                ok = 0;
          /* That was the membership test. ok true means we have found a new
           * rep */
          if (ok) {
            if (pthere) {
              ad = fpt[bno];
              while (*tailad[ad] == *trad[ad])
                ad = fpt[ad];
              tree[trptr] = *tailad[ad];
            }
            else {
              pthere = 1;
              ad = bno;
              tree[trptr] = pt;
              if (expg) {
                expp[pt] = p1;
                p1 = pst;
                p2 = p1 + npt;
                pst += (2 * npt);
                nexp++;
                if (pst > pend) {
                  fprintf(stderr, "Out of perm space. Increase PSP.\n");
                  return (-1);
                }
              }
            }
            *(trad[ad] + 1) = trptr;
            trad[ad] = tree + trptr;
            trptr += 2;
            for (i = fpt[ad]; i != -1; i = fpt[i]) {
              *(trad[i] + 1) = -1;
              tree[trptr] = *tailad[i];
              trad[i] = tree + trptr;
              trptr += 2;
            }
            cind++;
            tree[trptr] = cind;
            trptr++;
            if (trptr > trsp) {
              fprintf(stderr, "Out of tree space. Increase TRSP.\n");
              return (-1);
            }
            if (dcr)
              for (i = fnt; i != -1; i = fpt[i])
                fprintf(opx, "%5d", *trad[i]);
            /* if cr, then we output each coset rep */
            if (cr) {
              if (npt >= 1000) {
                for (i = 1; i <= npt; i++)
                  fprintf(op, "%5d", image(i));
                fprintf(op, "\n");
              }
              else {
                for (i = 1; i <= npt; i++)
                  fprintf(op, "%4d", image(i));
                fprintf(op, "\n");
              }
            }
            if (cind == ind)
              goto nextbno;
          } /* if (ok) */
          coind++;
        } /* while (coind<=oind) */
        (*gorb)++;
        gorb[*gorb] = pt;
      } /* if (svpt!=0 ... */
    }   /* for (pt=1;pt<=npt;... */
    fprintf(stderr, "Premature end of search.\n");
    return (-1);
  nextbno:
    if (expg)
      for (i = 1; i <= npt; i++)
        svgptr[bno][i] = expp[i];
    sad -= 2;
    oind = ind;
  } /* main bno loop */
  if (dcr)
    fclose(opx);
  if (cr)
    fclose(op);
  for (i = fnt; i != -1; i = fpt[i])
    *(trad[i] + 1) = -1;
  printf("Tree space used = %d.\n", trptr);
  if (expg)
    printf("Number of perms expanded = %d.\n", nexp);
  return (0);
}

void advance(void)
/* This advances the element h in the search through elements gh */
{
  int ad, k, *p, *q;
  ad = lnt;
  while ((k = *(tailad[ad] + 1)) == -1)
    ad = bpt[ad];
  q = tree + k;
  ucp = cpad[ad] - 1;
  while (ad != -1) {
    tailad[ad] = q;
    cpad[ad] = ucp + 1;
    if (expg) {
      p = svgptr[ad][*q];
      if (p != stop) {
        ucp++;
        *ucp = p;
      }
    }
    else
      addsvf(*q, svgptr[ad]);
    q += 2;
    ad = fpt[ad];
  }
}
