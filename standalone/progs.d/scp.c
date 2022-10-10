#include "defs.h"
#include "permfns.h"

extern char mult, subgp, sgc, sgstr[], inf0[], inf1[], inf2[], inf3[], inf4[],
    outf[], outft[];
extern short perm[], sv[], cp[], fpt[], orb[], intpow[], base[], lorb[],
    pno[], *pptr[], *svptr[], intno[], ngno[], power[], wt[], co[], rel[],
    d1[], d2[], igno[], pinv[], rel[], pwt[], *sv2ptr[], crep[], crepinv[],
    dcrep[], dcrepinv[], tp[], mp, mb, mpt, mexp;
extern int psp, svsp;
char       norm;
short npt, npt1, nb, exp, prime, hgen, coeff, intexp, nwt, stig, *itp, ngads,
    mxp;
FILE *ip, *ipcr, *ipkp, *op, *opy;

int scprog(void)
{
  short i, j, k, l, m, n, stconj, pm1, ncr, crct, ndc, dcct, lo, intexsk,
      stpt, lpt, stint, olen, pt, ad, *p1, *ip1, *p2, *ip2, *pint, *ptr, olo,
      **vsvptr;
  int  quot;
  char ingp, igth;

  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &exp, &nb, &l);
  if (npt > mpt) {
    fprintf(stderr, "npt too big. Increase MPT.\n");
    return (-1);
  }
  if (nb >= mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (nb * npt * 2 > svsp) {
    fprintf(stderr, "svsp too small. Increase SVSP.\n");
    return (-1);
  }
  if (exp >= mexp) {
    fprintf(stderr, "exp too big. Increase MEXP.\n");
    return (-1);
  }
  if (l != 2) {
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
  if (2 * exp > mxp) {
    fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + i * npt1 - 1;
  for (i = 1; i <= nb; i++)
    svptr[i] = sv + (i - 1) * npt - 1;
  readpsv(0, nb, exp, svptr);
  for (i = exp; i >= 1; i--)
    fscanf(ip, "%hd", pwt + i);
  fscanf(ip, "%hd%hd", &prime, &ngads);
  if (2 * exp + ngads > mxp) {
    fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  for (i = 1; i <= exp; i++)
    fscanf(ip, "%hd", power + i);
  for (i = 1; i <= exp; i++) {
    j = 2 * (i - 1);
    ngno[i] = j;
    igno[j + 1] = i;
  }
  k = 2 * exp - 1;
  for (i = 1; i <= ngads; i++) {
    l = i + k;
    readvec(pptr[l], 1);
    m = pptr[l][npt1];
    ngno[m] = l;
  }
  fclose(ip);

  stconj = 2 * exp + ngads - 1;
  itp = tp + npt1;
  if (4 * exp + ngads > mxp) {
    fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  if (subgp) {
    stig = stconj + 2 * exp;
    for (i = 1; i <= nb; i++)
      sv2ptr[i] = sv + (nb + i - 1) * npt - 1;
    if (subgp > 1) {
      strcpy(inf3, inf0);
      strcat(inf3, "sg");
      sgstr[0] = sgc;
      strcat(inf3, sgstr);
    }
    if (rsgp() == -1)
      return (-1);
    igth = 1;
  }
  else
    vsvptr = svptr;

  if ((ip = fopen(inf2, "r")) == 0) {
    printf("Cannot open %s.\n", inf2);
    return (-1);
  }
  fscanf(ip, "%hd%hd", &i, &ndc);
  while (getc(ip) != '\n')
    ;
  if (i != npt) {
    printf("dcr has npt false.\n");
    return (-1);
  }

  /* If subgp is true, then we will have to update the dcr file, so we print
     it out on to a temporary file as we go along.
  */
  if (subgp) {
    opy = fopen(outft, "w");
    fprintf(opy, "%4d%4d%4d%4d\n", npt, ndc, 0, 0);
  }
  op = fopen(outf, "w");
  pm1 = prime - 1;
  setpinv();

  for (dcct = 1; dcct <= ndc; dcct++) {
    if (subgp > 1 && igth == 0) {
      strcpy(inf3, inf0);
      strcat(inf3, "sg");
      sgstr[0] = sgc;
      strcat(inf3, sgstr);
      ipkp = ip;
      if (rsgp() == -1)
        return (-1);
      ip = ipkp;
      vsvptr = sv2ptr;
      igth = 1;
    }
    else if (subgp)
      vsvptr = sv2ptr;
    readvec(dcrep, 0);
    invert(dcrep, dcrepinv);
    fscanf(ip, "%hd", &lo);
    intexp = 0;
    olo = lo;
    while (lo % prime == 0) {
      intexp++;
      lo /= prime;
    }
    if (lo != 1)
      printf("Warning. lo was not a power of p.\n");
    intexp = exp - intexp;
    printf("dcct,intexp=%d,%d.\n", dcct, intexp);
    /* intexp is the exponent of Q = P ^ gPg(-1), for the current dcrep g */
    if (intexp == 0 || (intexp == 1 && mult)) {
      if (subgp) {
        if (npt >= 1000)
          for (n = 1; n <= npt; n++)
            fprintf(opy, "%5d", dcrep[n]);
        else
          for (n = 1; n <= npt; n++)
            fprintf(opy, "%4d", dcrep[n]);
        fprintf(opy, "   %4d\n", olo);
      }
      if (mult == 0)
        fprintf(op, "%4d\n", 0);
      continue;
    }
    /* Now we compute the gens g(-1)h(i)g, where h(i) are the PCP gens of P */
    for (i = 1; i <= exp; i++) {
      p1 = pptr[ngno[i]];
      p2 = pptr[stconj + i];
      for (n = 1; n <= npt; n++)
        p2[n] = dcrep[p1[dcrepinv[n]]];
    }
    stint = stconj + exp;
    pint = pptr[stint];
    intexsk = 0;
    for (i = 1; i <= exp; i++)
      intno[i] = 0;
    lpt = 1;
    stpt = 1;
    pm1 = prime - 1;
    /* Now we search through the elements of g(-1)Pg, testing for membership
       of P, until we have found the intexp gens of g(-1)Qg.
    */
    for (i = 1; i <= exp; i++) {
      olen = 1;
      cp[1] = stconj + i;
      for (j = 1; j <= exp; j++)
        co[j] = 0;
      while (1) {
        *cp = olen;
        ingp = 1;
        for (k = 1; k <= nb; k++) {
          pt = image(base[k]);
          ptr = vsvptr[k];
          if (ptr[pt] == 0) {
            ingp = 0;
            ad = stpt;
            break;
          }
          if (k < nb)
            addsv(pt, ptr);
        }
        if (ingp)
          break;
        while (co[ad] == pm1) {
          olen -= pm1;
          co[ad] = 0;
          ad = fpt[ad];
        }
        if (ad == i) {
          lpt = i;
          fpt[i] = i + 1;
          break;
        }
        olen++;
        co[ad]++;
        cp[olen] = stconj + ad;
      }
      if (ingp) {
        intexsk++;
        pint += npt1;
        *cp = olen;
        for (n = 1; n <= npt; n++)
          pint[n] = image(n);
        pint[npt1] = intexsk;
        intno[i] = stconj + 2 * intexsk - 1;
        if (intexsk == intexp)
          break;
        if (stpt == i) {
          lpt = i + 1;
          stpt = lpt;
        }
        else
          fpt[lpt] = i + 1;
      }
    }
    /* Search is complete */
    if (intexsk < intexp) {
      fprintf(stderr, "Intersection error.\n");
      return (-1);
    }
    printf("Found intersection.\n");

    /* If subgp, then we will now modify g, until H ^ gHg- contains Q */
    if (subgp)
      for (m = subgp; m >= 1; m--) {
        if (m == subgp)
          for (i = 1; i <= npt; i++) {
            tp[i] = i;
            itp[i] = i;
          }
        if (subgp > 1) {
          strcpy(inf4, inf0);
          strcat(inf4, "cr");
          strcat(inf4, sgstr);
          if (m > 1) {
            sgstr[0]--;
            strcpy(inf3, inf0);
            strcat(inf3, "sg");
            strcat(inf3, sgstr);
            ipkp = ip;
            if (rsgp() == -1)
              return (-1);
            ip = ipkp;
            igth = 0;
            vsvptr = sv2ptr;
          }
          else
            vsvptr = svptr;
        }
        else
          vsvptr = svptr;
        if ((ipcr = fopen(inf4, "r")) == 0) {
          fprintf(stderr, "Cannot open %s.\n", inf4);
          return (-1);
        }
        ingp = 0;
        fscanf(ipcr, "%hd%hd", &n, &ncr);
        while (getc(ipcr) != '\n')
          ;
        if (n != npt) {
          fprintf(stderr, "inf4 has npt wrong.\n");
          return (-1);
        }
        for (crct = 0; crct <= ncr; crct++) {
          if (crct == 0)
            for (n = 1; n <= npt; n++)
              crep[n] = n;
          else
            for (n = 1; n <= npt; n++)
              fscanf(ipcr, "%hd", crep + n);
          invert(crep, crepinv);
          for (i = 1; i <= intexp; i++) {
            p1 = pptr[stint + i];
            *cp = 0;
            for (j = 1; j <= nb; j++) {
              pt = image(crep[tp[p1[itp[crepinv[base[j]]]]]]);
              ptr = vsvptr[j];
              if (ptr[pt] == 0)
                goto nextcr;
              if (j < nb)
                addsv(pt, ptr);
            }
          }
          ingp = 1;
          printf("Got intersection in subgp %d. crct=%d.\n", m - 1, crct);
          for (n = 1; n <= npt; n++)
            tp[n] = crep[tp[n]];
          invert(tp, itp);
          fclose(ipcr);
          break;
        nextcr:;
        }
        if (ingp == 0) {
          fprintf(stderr, "Cannot get intersection in subgp %d.\n", m - 1);
          return (-1);
        }
      }

    /* Now we reconjugate the gens of g(-1)Qg by g(-1) to get gens of Q */
    for (i = 1; i <= intexp; i++) {
      p1 = pptr[stconj + 2 * i - 1];
      p2 = pptr[stint + i];
      for (n = 1; n <= npt; n++)
        p1[n] = dcrepinv[p2[dcrep[n]]];
      p1[npt1] = p2[npt1];
      invert(p1, p1 + npt1);
    }
    /* If subgp, then update dcr */
    if (subgp) {
      if (npt >= 1000)
        for (n = 1; n <= npt; n++) {
          dcrep[n] = tp[dcrep[n]];
          fprintf(opy, "%5d", dcrep[n]);
        }
      else
        for (n = 1; n <= npt; n++) {
          dcrep[n] = tp[dcrep[n]];
          fprintf(opy, "%4d", dcrep[n]);
        }
      fprintf(opy, "   %4d\n", olo);
      invert(dcrep, dcrepinv);
    }
    norm = intexp == exp;
    /* If norm, then we do not need to compute the PCP for Q */
    if (norm) {
      fprintf(op, "%4d\n", intexp);
      if (mult == 0) {
        for (i = exp; i >= 1; i--) {
          wt[i] = pwt[i];
          fprintf(op, "%4d", wt[i]);
        }
        fprintf(op, "\n");
      }
      goto outconj;
    }

    /* Now we compute the PCP for Q. This is similar to the algorithm in pcrun
     */
    for (i = 1; i <= intexp; i++) {
      wt[i] = 1;
      d1[i] = 0;
      d2[i] = 0;
    }
  restart:
    if (mult)
      nwt = 2;
    for (i = intexp; i >= 2; i--) {
      p1 = pptr[stconj + 2 * i - 1];
      ip1 = p1 + npt1;
      for (j = intexp; j > i; j--) {
        if (mult == 0)
          nwt = wt[i] + wt[j];
        p2 = pptr[stconj + 2 * j - 1];
        ip2 = p2 + npt1;
        for (n = 1; n <= npt; n++) {
          pt = p2[p1[ip2[ip1[n]]]];
          tp[n] = pt;
          itp[pt] = n;
        }
        if ((n = expint(intexp + 1 - i, intexp + 1 - j, nwt)) > 0)
          goto restart;
        if (n == -1)
          return (-1);
      }
    }
    if (mult == 0)
      for (i = intexp; i >= 2; i--) {
        nwt = wt[i] + 1;
        p1 = pptr[stconj + 2 * i - 1];
        ip1 = p1 + npt1;
        for (n = 1; n <= npt; n++) {
          pt = n;
          for (m = 1; m <= prime; m++)
            pt = p1[pt];
          tp[n] = pt;
          itp[pt] = n;
        }
        if ((n = expint(intexp + 1 - i, intexp + 1 - i, nwt)) > 0)
          goto restart;
        if (n == -1)
          return (-1);
      }
    fprintf(op, "%4d\n", intexp);
    if (mult == 0) {
      for (i = intexp; i >= 1; i--)
        fprintf(op, "%3d", wt[i]);
      fprintf(op, "\n");
    }
    /* We output the PCP gens of Q followed by those of g(-1)Qg */
    for (i = intexp; i >= 1; i--) {
      p1 = pptr[stconj + 2 * i - 1];
      p2 = tp;
      ptr = p1 + 2 * npt1;
      while (p1 < ptr)
        *(++p2) = *(++p1);
      express(tp, rel, 0);
      l = *rel;
      for (n = 0; n <= l; n++)
        fprintf(op, "%4d", rel[n]);
      fprintf(op, "\n");
    }
  outconj:
    for (i = intexp; i >= 1; i--) {
      p1 = pptr[stconj + 2 * i - 1];
      for (n = 1; n <= npt; n++) {
        pt = dcrep[p1[dcrepinv[n]]];
        tp[n] = pt;
        itp[pt] = n;
      }
      express(tp, rel, 0);
      l = *rel;
      for (n = 0; n <= l; n++)
        fprintf(op, "%4d", rel[n]);
      fprintf(op, "\n");
    }
    if (norm)
      continue;
    if (mult == 0) {
      for (i = intexp; i >= 1; i--)
        fprintf(op, "%3d", d1[i]);
      fprintf(op, "\n");
      for (i = intexp; i >= 1; i--)
        fprintf(op, "%3d", d2[i]);
      fprintf(op, "\n");
    }
    for (i = intexp; i >= 2; i--) {
      p1 = pptr[stconj + 2 * i - 1];
      ip1 = p1 + npt1;
      for (j = intexp; j > i; j--) {
        p2 = pptr[stconj + 2 * j - 1];
        ip2 = p2 + npt1;
        for (n = 1; n <= npt; n++) {
          pt = p2[p1[ip2[ip1[n]]]];
          tp[n] = pt;
          itp[pt] = n;
        }
        expint(intexp + 1 - i, intexp + 1 - j, 0);
      }
    }
    if (mult == 0)
      for (i = intexp; i >= 2; i--) {
        p1 = pptr[stconj + 2 * i - 1];
        ip1 = p1 + npt1;
        for (n = 1; n <= npt; n++) {
          pt = n;
          for (m = 1; m <= prime; m++)
            pt = p1[pt];
          tp[n] = pt;
          itp[pt] = n;
        }
        expint(intexp + 1 - i, intexp + 1 - i, 0);
      }
  }
  fprintf(op, "%d\n", -1);
  if (subgp) {
    fclose(opy);
    fclose(op);
    fclose(ip);
    ip = fopen(outft, "r");
    op = fopen(inf2, "w");
    while ((i = getc(ip)) != -1)
      putc(i, op);
    fclose(ip);
    unlink(outft);
  }
  return (0);
}

int rsgp(void)
/* Used to read in intermediate subgroups when subgp is true */
{
  short i, np, k, l, n, *ptr;
  if ((ip = fopen(inf3, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf3);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &i, &np, &k, &l);
  if (i != npt || k != nb || l <= 0) {
    fprintf(stderr, "Intermediate group has wrong input format.\n");
    return (-1);
  }
  seeknln();
  seeknln();
  seeknln();
  if (4 * exp + 2 * np + ngads > mxp) {
    fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  for (i = 1; i <= np; i++) {
    readvec(pptr[stig + i + 1], 1);
    invert(pptr[stig + i + 1], pptr[stig + i]);
  }
  for (i = 1; i <= nb; i++) {
    readvec(sv2ptr[i], 0);
    for (n = 1; n <= npt; n++) {
      ptr = sv2ptr[i] + n;
      if (*ptr > 0)
        *ptr = stig + (*ptr + 1) / 2;
    }
  }
  fclose(ip);
  return (0);
}

int expint(int i, int j, int nwt)
/* This plays the same role for Q as express in pcscfns.c plays for P */
{
  short l, m, n, ino, pt, p, *g, *ig;
  char  diff;
  l = 1;
  while (1) {
    firstgen(tp, &hgen, &coeff);
    if (hgen == 0)
      break;
    g = pptr[intno[hgen]];
    ig = g + npt1;
    if (g == 0) {
      printf("Subint error.\n");
      return (-1);
    }
    ino = g[npt1];
    if ((wt[ino] < nwt) || (wt[ino] == nwt && d1[ino] == 0)) {
      wt[ino] = nwt;
      intpow[ino] = pinv[coeff];
      diff = 0;
      for (n = 1; n <= npt; n++)
        if (g[n] != tp[n]) {
          diff = 1;
          break;
        }
      if (diff) {
        for (n = 1; n <= npt; n++) {
          pt = tp[n];
          g[n] = pt;
          ig[pt] = n;
        }
        for (n = 1; n <= intexp; n++) {
          d1[n] = 0;
          d2[n] = 0;
        }
        return (ino);
      }
      else {
        d1[ino] = i;
        d2[ino] = j;
      }
    }
    l++;
    rel[l] = 1 + intexp - ino;
    l++;
    p = (coeff * intpow[ino]) % prime;
    rel[l] = p;
    for (n = 1; n <= npt; n++) {
      pt = itp[n];
      for (m = 1; m <= p; m++)
        pt = g[pt];
      tp[pt] = n;
    }
    invert(tp, itp);
  }
  if (nwt)
    return (0);
  m = 1 + intexp - rel[l - 1];
  if (d1[m] == i && d2[m] == j)
    rel[0] = rel[l - 1];
  else
    rel[0] = 0;
  rel[1] = l - 1;
  for (n = 0; n <= l; n++)
    fprintf(op, "%4d", rel[n]);
  fprintf(op, "\n");
  return (0);
}
