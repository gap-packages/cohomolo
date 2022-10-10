#include "defs.h"
#include "permfns.h"

extern char mult, conv, gens, inf1[], inf2[], outf1[], outf2[], outf3[],
    outf4[];
extern short perm[], sv[], cp[], fpt[], orb[], igno[], base[], lorb[], pno[],
    *pptr[], *svptr[], gno[], ngno[], power[], wt[], d1[], d2[], facord[],
    pinv[], rel[], *relptr[], mp, mpt, mb, mexp;
extern int psp, svsp;
short      npt, np, npt1, nb, exp, prime;
FILE *     ip, *op;

int resetsv(void)
/* This recomputes the Schreier vectors, after a change of generators */
{
  short i, j, lo, bno, b;
  for (bno = 1; bno <= nb; bno++)
    if (lorb[bno] > 1) {
      *pno = 0;
      lo = 0;
      b = base[bno];
      for (i = 1; i <= np; i++) {
        j = gno[i];
        if (pptr[j][npt1] >= bno) {
          (*pno)++;
          pno[*pno] = j;
          lo = orbitsv(b, svptr[bno], lo);
        }
      }
      if (lo != lorb[bno]) {
        fprintf(stderr, "Orbit length wrong.\n");
        return (-1);
      }
    }
  return (0);
}

void setfixb(short * p)
{
  short bno;
  bno = 1;
  while (bno <= nb && p[base[bno]] == base[bno])
    bno++;
  p[npt1] = bno;
}

int pcprog(void)
{
  short i, j, k, l, m, n, mxp, ngads, class, hgen, coeff, nf, nf1, nf2, pt,
      *pk, *pnf, *pnf1, *pnf2, *ipnf, *ipnf1, *p1, *p2, *ip1, *ip2, nwt;
  int  quot;
  char diff;
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &np, &nb, &l);
  if (npt > mpt) {
    fprintf(stderr, "npt too big. Increase MPT.\n");
    return (-1);
  }
  if (nb >= mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (nb * npt > svsp) {
    fprintf(stderr, "Out of sv space. Increase SVSP.\n");
    return (-1);
  }
  if (l != 4) {
    fprintf(stderr, "Wrong input format.\n");
    return (-1);
  }
  readbaselo(nb, base, lorb);
  npt1 = npt + 1;
  quot = psp / npt1;
  if (quot > mp)
    quot = mp;
  mxp = quot;
  mxp = 2 * (mxp / 2);
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + i * npt1 - 1;
  for (i = 1; i <= nb; i++)
    svptr[i] = sv + (i - 1) * npt - 1;
  if (2 * (np + 1) > mxp) {
    fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  readpsv(0, nb, np, svptr);
  fscanf(ip, "%hd", &prime);
  for (i = 1; i <= np; i++) {
    fscanf(ip, "%hd", facord + i);
    gno[i] = 2 * (i - 1);
    igno[2 * i - 1] = i;
  }
  fclose(ip);

  /* Since we will be changing generators, we link the perm nos with fpt.
     gno[i] will be the current gen no i in the subnormal series for P.
  */
  nf = 2 * np;
  for (i = nf; i < mxp; i += 2)
    fpt[i] = i + 2;
  fpt[mxp - 2] = -1;

  /* Now we start to convert the subnormal series to a central series */
  for (i = 1; i <= np; i++) {
    p1 = pptr[gno[i]];
    m = facord[i];
    /* If m>prime, we introduce p-th powers of this generator, until all
       indices in the series have order p.
    */
    while (m > prime) {
      pnf = pptr[nf];
      ipnf = pnf + npt1;
      for (n = 1; n <= npt; n++) {
        pt = n;
        for (l = 1; l <= prime; l++)
          pt = p1[pt];
        pnf[n] = pt;
        ipnf[pt] = n;
      }
      for (n = np + 1; n >= i + 1; n--) {
        l = gno[n - 1];
        gno[n] = l;
        igno[l + 1] = n;
        facord[n + 1] = facord[n];
      }
      setfixb(pnf);
      facord[i + 1] = prime;
      m /= prime;
      facord[i] = m;
      gno[i] = nf;
      igno[nf + 1] = i;
      np++;
      p1 = pnf;
      nf = fpt[nf];
      if (nf == -1) {
        fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
        return (-1);
      }
      resetsv();
    }
    ip1 = p1 + npt1;
    /* In the main bit of the algorithm, we compute commutator [i,j], and
       change j until this lies in G(j-1).
    */
    for (j = 1; j <= (i - 2); j++) {
      p2 = pptr[gno[j]];
      while (1) {
        ip2 = p2 + npt1;
        pnf = pptr[nf];
        ipnf = pnf + npt1;
        for (n = 1; n <= npt; n++) {
          l = p1[p2[ip1[ip2[n]]]];
          pnf[n] = l;
          ipnf[l] = n;
        }
        setfixb(pnf);
        firstgen(pnf, &hgen, &coeff);
        k = hgen;
        if (k >= i) {
          fprintf(stderr, "Series is not subnormal.\n");
          return (-1);
        }
        if (k < j)
          break;
        nf1 = fpt[nf];
        pnf1 = pptr[nf1];
        ipnf1 = pnf1 + npt1;
        nf2 = fpt[nf1];
        if (nf1 == -1 || nf2 == -1) {
          fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
          return (-1);
        }
        pk = pptr[gno[k]];
        for (n = 1; n <= npt; n++) {
          pt = n;
          for (l = 1; l <= coeff; l++)
            pt = pk[pt];
          pnf1[n] = pt;
          ipnf1[pt] = n;
        }
        setfixb(pnf1);
        fpt[gno[k]] = nf2;
        nf2 = gno[k];
        fpt[nf1] = nf2;
        pnf2 = pptr[nf2];
        while (1) {
          for (n = 1; n <= npt; n++)
            pnf2[n] = pnf[ipnf1[n]];
          firstgen(pnf2, &hgen, &coeff);
          if (hgen >= k) {
            fprintf(stderr, "Impossible error!\n");
            return (-1);
          }
          if (hgen < j) {
            for (m = k; m > j; m--) {
              l = gno[m - 1];
              gno[m] = l;
              igno[l + 1] = m;
            }
            gno[j] = nf1;
            p2 = pnf1;
            igno[nf1 + 1] = j;
            fpt[nf] = nf2;
            resetsv();
            break;
          }
          for (m = k; m >= (hgen + 2); m--) {
            l = gno[m - 1];
            gno[m] = l;
            igno[l + 1] = m;
          }
          if (pptr[gno[hgen]][npt1] >= pnf1[npt1]) {
            l = gno[hgen];
            gno[hgen + 1] = l;
            igno[l + 1] = hgen + 1;
          }
          else {
            gno[hgen + 1] = nf1;
            igno[nf1 + 1] = hgen + 1;
            nf1 = gno[hgen];
            fpt[nf1] = nf2;
            fpt[nf] = nf1;
            pnf1 = pptr[nf1];
          }
          pk = pptr[gno[hgen] + 1];
          for (n = 1; n <= npt; n++) {
            pt = n;
            for (l = 1; l <= coeff; l++)
              pt = pk[pt];
            pt = ipnf1[pt];
            pnf1[pt] = n;
          }
          ipnf1 = pnf1 + npt1;
          invert(pnf1, ipnf1);
          setfixb(pnf1);
          if (hgen == j) {
            p2 = pnf1;
            gno[j] = nf1;
            igno[nf1 + 1] = j;
            fpt[nf] = nf2;
            resetsv();
            break;
          }
          k = hgen;
        } /* while(1) */
      }   /* while(1)  */
    }     /* j loop  */
  }       /* i loop */
  printf("Central series found.\n");

  /* We are now ready to compute the pcp, with definitions (d1,d2) and weigths
     (wt). This involves two stages. In the first, we compute all commutators
     and powers, to determine what the weights and definitios are, and change
     some of the generators if necessary. We do not change the Schreier vector
     generators - the new generators are used as PCP generators. Whenever we
     change a generator, we have to completely restart this process. In the
     second stage, the PCP is computed and output.
  */
  exp = np;
  setpinv();
  if (exp >= mexp) {
    fprintf(stderr, "exp too big. Increase MEXP.\n");
    return (-1);
  }
  for (i = 1; i <= exp; i++) {
    ngno[i] = gno[i];
    wt[i] = 1;
    power[i] = 1;
    d1[i] = 0;
    d2[i] = 0;
  }
  ngads = 0;
  class = 1;

restart:
  for (i = exp; i > 1; i--)
    for (j = exp; j >= i; j--) {
      pnf = pptr[nf];
      ipnf = pnf + npt1;
      p1 = pptr[ngno[i]];
      ip1 = p1 + npt1;
      if (i == j) {
        if (mult)
          continue;
        nwt = wt[i] + 1;
        for (n = 1; n <= npt; n++) {
          pt = n;
          for (l = 1; l <= prime; l++)
            pt = p1[pt];
          pnf[n] = pt;
          ipnf[pt] = n;
        }
      }
      else {
        nwt = wt[i] + wt[j];
        p2 = pptr[ngno[j]];
        ip2 = p2 + npt1;
        for (n = 1; n <= npt; n++) {
          pt = p2[p1[ip2[ip1[n]]]];
          pnf[n] = pt;
          ipnf[pt] = n;
        }
      }
      hgen = express(pnf, rel, nwt);
      if (hgen > 0) {
        wt[hgen] = nwt;
        if (nwt > class)
          class = nwt;
        l = ngno[hgen];
        diff = 0;
        p1 = pptr[l];
        for (k = 1; k <= npt; k++)
          if (pnf[k] != p1[k]) {
            diff = 1;
            break;
          }
        if (diff) {
          if (l != gno[hgen]) {
            fpt[l] = fpt[nf];
            fpt[nf] = l;
          }
          else
            ngads++;
          ngno[hgen] = nf;
          pptr[nf][npt1] = hgen;
          nf = fpt[nf];
          if (nf == -1) {
            fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
            return (-1);
          }
          for (k = 1; k <= exp; k++) {
            d1[k] = 0;
            d2[k] = 0;
          }
          goto restart;
        }
        else {
          d1[hgen] = exp + 1 - i;
          d2[hgen] = exp + 1 - j;
        }
      }
    }

  /* Now we are ready to compute the PCP. We output it as we go along. */
  op = fopen(outf1, "w");
  fprintf(op, "%4d%4d%4d%4d\n", npt, np, nb, 2);
  printbaselo(nb, base, lorb);
  *pno = 0;
  *gno = np;
  printpsv(nb, gno, svptr);
  for (i = exp; i > 0; i--)
    fprintf(op, "%4d", wt[i]);
  fprintf(op, "\n");
  fprintf(op, "%4d%4d\n", prime, ngads);
  for (i = 1; i <= exp; i++)
    fprintf(op, "%4d", power[i]);
  fprintf(op, "\n");
  for (i = 1; i <= exp; i++) {
    j = ngno[i];
    if (j != gno[i])
      printvec(pptr[j], 1);
  }
  fclose(op);
  /* Output pcp gens as permutations if required. */
  if (gens) {
    op = fopen(outf4, "w");
    fprintf(op, "%4d%4d%4d%4d\n", npt, exp, 0, 0);
    for (i = exp; i >= 1; i--) {
      j = ngno[i];
      printvec(pptr[j], 0);
    }
    fclose(op);
  }
  op = fopen(outf2, "w");
  fprintf(op, "%4d%4d%4d%4d%4d%4d\n", prime, exp, exp, exp - 1, class, mult);
  for (i = exp; i > 0; i--)
    fprintf(op, "%4d", wt[i]);
  fprintf(op, "\n");
  for (i = exp; i > 0; i--)
    fprintf(op, "%4d", d1[i]);
  fprintf(op, "\n");
  for (i = exp; i > 0; i--)
    fprintf(op, "%4d", d2[i]);
  fprintf(op, "\n");
  for (i = 2; i < exp; i++)
    for (j = 1; j < i; j++) {
      pnf = pptr[nf];
      ipnf = pnf + npt1;
      p1 = pptr[ngno[exp + 1 - i]];
      ip1 = p1 + npt1;
      p2 = pptr[ngno[exp + 1 - j]];
      ip2 = p2 + npt1;
      for (n = 1; n <= npt; n++) {
        pt = p2[p1[ip2[ip1[n]]]];
        pnf[n] = pt;
        ipnf[pt] = n;
      }
      express(pnf, rel + 1, 0);
      l = *(rel + *(rel + 1));
      m = 1 + exp - l;
      if (d1[m] == i && d2[m] == j)
        *rel = l;
      else
        *rel = 0;
      for (k = 0; k <= *(rel + 1) + 1; k++)
        fprintf(op, "%4d", rel[k]);
      fprintf(op, "\n");
    }
  for (i = 1; i < exp; i++) {
    pnf = pptr[nf];
    ipnf = pnf + npt1;
    p1 = pptr[ngno[exp + 1 - i]];
    ip1 = p1 + npt1;
    for (n = 1; n <= npt; n++) {
      pt = n;
      for (l = 1; l <= prime; l++)
        pt = p1[pt];
      pnf[n] = pt;
      ipnf[pt] = n;
    }
    express(pnf, rel + 1, 0);
    l = *(rel + *(rel + 1));
    m = 1 + exp - l;
    if (d1[m] == i && d2[m] == i)
      *rel = l;
    else
      *rel = 0;
    for (k = 0; k <= *(rel + 1) + 1; k++)
      fprintf(op, "%4d", rel[k]);
    fprintf(op, "\n");
  }
  fclose(op);
  if (conv) {
    if ((ip = fopen(inf2, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf2);
      return (-1);
    }
    op = fopen(outf3, "w");
    fscanf(ip, "%hd%hd%hd%hd", &k, &l, &m, &n);
    if (k != npt || n <= 0) {
      fprintf(stderr, "inf2 has invalid npt or format.\n");
      return (-1);
    }
    for (i = 1; i <= 3; i++)
      while (getc(ip) != '\n')
        ;
    fprintf(op, "%3d\n", l);
    pnf = pptr[nf];
    ipnf = pnf + npt1;
    for (i = 1; i <= l; i++) {
      readvec(pnf, 1);
      invert(pnf, ipnf);
      express(pnf, rel, 0);
      fprintf(op, "%4d  ", *rel);
      for (j = 1; j <= *rel; j++)
        fprintf(op, "%4d", rel[j]);
      fprintf(op, "\n");
    }
  }
  return (0);
}
