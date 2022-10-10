#include <stdio.h>
#define RD(A, B) (done[(B) / 32][(A)] & (one << ((B) % 32)))
#define WD(A, B) done[(B) / 32][(A)] |= (one << ((B) % 32))
extern char inf[], outf[], outfg[], gap, firstnew;
extern int  imsp[], *done[], *imcos[];
extern int  perm[], sv[], cp[], orb[], base[], lorb[], pno[], *pptr[],
    *svptr[], inv[], rno[], mp, mrel, dnwds, mpt, mb;
extern int psp, svsp, space;
int nr, endr, maxcos, *fpt, *bpt, *def, scanno, nscan, nelim, bno, lo, ccos,
    lastd, cind, nfree, stcr, endcr, fcos, bcos, lcl, np2, *rel, one = 1;
int   npt;
int   rsp;
char  fullsc, clsd, lkah;
FILE *fopen(), *ip, *op;
/* Many of the variables are as in program tcrun */

int grprog(void)
{
  int   i, j, k, l, m, n, b, l1, l2, nb, rn1, nperms, npt1, mxp, min, mini;
  int   spn, quot;
  float a;
  ip = fopen(inf, "r");
  if (ip == 0) {
    fprintf(stderr, "Cannot open %s\n", inf);
    return (-1);
  }
  fscanf(ip, "%d%d%d%d", &npt, &nperms, &nb, &k);
  if (npt > mpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    return (-1);
  }
  if (nb >= mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (nb * npt > svsp) {
    fprintf(stderr, "Out of sv space.Increase SVSP.\n");
    return (-1);
  }
  if (k <= 0) {
    fprintf(stderr, "Wrong input format.\n");
    return (-1);
  }
  npt1 = npt + 1;
  np2 = 2 * nperms - 1;
  quot = psp / npt1;
  if (quot > mp)
    quot = mp;
  mxp = quot;
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + i * npt1 - 1;
  for (i = 1; i <= nb; i++)
    svptr[i] = sv + (i - 1) * npt - 1;
  if (np2 >= mxp) {
    fprintf(stderr, "nperms too big. Increase PSP (or MP).\n");
    return (-1);
  }
  readbaselo(nb, base, lorb);
  readpsv(0, nb, nperms, svptr);
  fclose(ip);
  for (i = 0; i < np2; i += 2) {
    inv[i] = i + 1;
    inv[i + 1] = i;
  }

  /* Group is read in. Now initialize coset table for enumaration. */
  printf("Input a,b,l1,l2.\nmaxcos = a*lorb+b\n");
  printf("Exit lookahead after scanning l1 or eliminating l2 cosets.\n");
  scanf("%f%d%d%d", &a, &b, &l1, &l2);
  mrel = dnwds * 32;
  maxcos = a * npt + b;
  spn = maxcos * (np2 + 4 + dnwds) + 1;
  if (spn > space) {
    fprintf(stderr, "Not enough space in coset table. Increase SPACE.\n");
    return (-1);
  }
  rel = imsp + spn;
  rsp = space - spn;
  for (i = 0; i <= np2; i++)
    imcos[i] = imsp + i * maxcos - 1;
  fpt = imcos[np2] + maxcos + 1;
  bpt = fpt + maxcos;
  def = bpt + maxcos;
  for (i = 0; i < dnwds; i++)
    done[i] = def + maxcos * (i + 1);

  if (gap) {
    op = fopen(outfg, "w");
    fprintf(op, "COHOMOLO.PermRels := [\n");
  }
  nr = 0;
  endr = -1;
  for (bno = nb; bno >= 1; bno--) {
    printf("bno=%d.\n", bno);
    /* We are going to do a coset enumeration of the group G[bno] (in the
       stabilizer chain), using genrators of G[bno+1] as subgroup. First we
       write in as many entries in the coset table as possible, making
       definitions of the cosets corresponding to the points in the orbit of
       bno under G[bno], using the Schreier vector sv[bno]. At the end of the
       enumeration, only these cosets should remain, and none of them should
       get eliminated.
    */
    lo = lorb[bno];
    rno[bno] = nr;
    maxcos = lo * a + b;
    rn1 = nr - 1;
    *pno = 0;
    for (i = 0; i <= np2; i++)
      for (j = 1; j <= maxcos; j++)
        imcos[i][j] = 0;
    for (i = 0; i < np2; i += 2)
      if (pptr[i][npt1] >= bno) {
        (*pno)++;
        pno[*pno] = i;
        if (pptr[i][npt1] > bno) {
          imcos[i][1] = 1;
          imcos[i + 1][1] = 1;
        }
      }
    if (orbitsv(base[bno], svptr[bno], 0) != lo) {
      fprintf(stderr, "Orbit length wrong!\n");
      return (-1);
    }
    if (lo == 1)
      continue;
    def[1] = -1;
    for (l = 1; l <= lo; l++) {
      n = orb[l];
      cp[n] = l;
      if ((k = svptr[bno][n]) != -1) {
        j = pptr[k][n];
        m = cp[j];
        imcos[k - 1][m] = l;
        imcos[k][l] = m;
        def[l] = k;
      }
    }
    /* Now we initialize the relator set, by putting in one relator for each
       generator in G[bno].
    */
    for (i = 0; i < np2; i += 2)
      if (pptr[i][npt1] == bno) {
        *cp = 1;
        cp[1] = i + 1;
        for (j = bno; j <= nb; j++)
          addsv(image(base[j]), svptr[j]);
        if (addrel() == -1)
          return (-1);
      }

    /* Since new relators will be added in the course of the enumeration, we
       want to avoid scanning the same coset under the same relator more than
       once. The array done records which pairs have been scanned in
       compressed form, using one bit for each coset-relator pair. WD(a,b)
       writes a one into the bit corresponding to coset a and relator b.
       RD(a,b) returns the value of this bit.
    */
    for (j = 0; j < dnwds; j++)
      for (i = 1; i <= lo; i++)
        done[j][i] = 0;
    for (i = 0; i <= rn1; i++)
      WD(1, i);
    /* Now we start the enumeration, as in tcrun  */
    lastd = lo;
    cind = lo;
    nfree = lo + 1;
    for (i = 1; i <= lo; i++)
      bpt[i] = i - 1;
    for (i = 0; i < maxcos; i++)
      fpt[i] = i + 1;
    fpt[lo] = 0;
    fpt[maxcos] = 0;
    while (1) {
      ccos = 1;
      lkah = 0;
      while (ccos != 0) {
        clsd = 1;
        endcr = -1;
        scanno = -1;
        while (endcr != endr) {
          scanno++;
          fullsc = 1;
          stcr = endcr + 2;
          endcr += (1 + rel[stcr - 1]);
          if (RD(ccos, scanno) == 0)
            scanrel();
          if (fullsc == 0) {
            clsd = 0;
            if (lkah == 0) {
              lcl = bpt[ccos];
              lkah = 1;
              nscan = 0;
              nelim = 0;
              /* printf("Entering lookahead.\n"); */
            }
          }
        }
        ccos = fpt[ccos];
        if (lkah) {
          nscan++;
          if (nscan > l1 || nelim >= l2 || ccos == 0) {
            if (cind !=
                maxcos) { /* printf("Exiting lookahead. cind=%d.\n",cind); */
              ccos = fpt[lcl];
              lkah = 0;
            }
            else
              ccos = 0;
          }
        }
      }
      if (cind == lo)
        break;
      /* The enumeration did not complete, so we add a new relator */
      if (firstnew)
        i = fpt[lo];
      else {
        min = 0;
        for (n = fpt[lo]; n != 0; n = fpt[n]) {
          *cp = 0;
          i = n;
          for (j = def[i]; j != -1; j = def[i]) {
            (*cp)++;
            cp[*cp] = j;
            i = imcos[j][i];
          }
          for (i = 1, j = *cp; i <= j; i++, j--) {
            if (i == j)
              cp[i] = inv[cp[i]];
            else {
              k = cp[i];
              cp[i] = inv[cp[j]];
              cp[j] = inv[k];
            }
          }
          for (i = bno; i <= nb; i++)
            addsv(image(base[i]), svptr[i]);
          if (min == 0 || *cp < min) {
            min = *cp;
            mini = n;
          }
        }
        i = mini;
      }
      *cp = 0;
      for (j = def[i]; j != -1; j = def[i]) {
        (*cp)++;
        cp[*cp] = j;
        i = imcos[j][i];
      }
      for (i = 1, j = *cp; i <= j; i++, j--) {
        if (i == j)
          cp[i] = inv[cp[i]];
        else {
          k = cp[i];
          cp[i] = inv[cp[j]];
          cp[j] = inv[k];
        }
      }
      for (i = bno; i <= nb; i++)
        addsv(image(base[i]), svptr[i]);
      if (addrel() == -1)
        return (-1);
    }
  }
  if (gap) {
    fprintf(op, "];\n");
    fclose(op);
  }
  op = fopen(outf, "w");
  fprintf(op, "%3d   ", nb);
  for (i = 1; i <= nb; i++)
    fprintf(op, "%4d", rno[i]);
  fprintf(op, "\n");
  m = -1;
  while (m != endr) {
    m++;
    l = rel[m];
    for (i = 0; i <= l; i++)
      fprintf(op, "%4d", rel[m + i]);
    fprintf(op, "\n");
    m += l;
  }
  return (0);
}

int addrel(void)
/* Adds a relator */
{
  int i, j, k;
  nr++;
  if (nr > mrel) {
    fprintf(stderr, "Too many relations. Increase MREL.\n");
    return (-1);
  }
  if (endr + *cp + 1 >= rsp) {
    fprintf(stderr, "Out of relation space. Increase SPACE.\n");
    return (-1);
  }
  endr++;
  rel[endr] = *cp;
  for (i = 1; i <= *cp; i++)
    rel[endr + i] = inv[cp[*cp + 1 - i]];
  rno[bno]++;
  endr += *cp;
  if (gap) {
    if (nr > 1)
      fprintf(op, ",\n");
    putc('[', op);
    for (i = 1; i <= *cp; i++) {
      j = cp[*cp + 1 - i];
      k = (j % 2 == 0) ? -(j / 2 + 1) : (j / 2 + 1);
      fprintf(op, "%d", k);
      if (i < *cp)
        putc(',', op);
      else
        putc(']', op);
    }
  }
  for (i = 1; i <= *cp; i++) {
    j = cp[*cp + 1 - i];
    k = (j % 2 == 0) ? -(j / 2 + 1) : (j / 2 + 1);
    printf("%2d", k);
    if (i < *cp)
      printf(" * ");
    else
      printf("\n");
  }
  return (0);
}

int scanrel(void)
{
  int  i, j, k, l, m;
  char comp;
  fcos = ccos;
  bcos = ccos;
  comp = 1;
  for (i = stcr; i <= endcr; i++) {
    k = imcos[rel[i]][fcos];
    if (k == 0) {
      comp = 0;
      break;
    }
    fcos = k;
  }
  if (comp) {
    if (fcos != bcos)
      coinc(fcos, bcos);
    WD(ccos, scanno);
    return (0);
  }
  stcr = i;
  for (i = endcr; i >= stcr; i--) {
    l = rel[i];
    m = inv[l];
    k = imcos[m][bcos];
    if (k == 0) {
      if (i == stcr) {
        k = imcos[l][fcos];
        if (k == 0) {
          imcos[l][fcos] = bcos;
          imcos[m][bcos] = fcos;
        }
        else if (k != bcos)
          coinc(k, bcos);
        WD(ccos, scanno);
        return (0);
      }
      if (lkah || nfree == 0) {
        fullsc = 0;
        return (0);
      }
      for (j = 0; j <= np2; j++)
        imcos[j][nfree] = 0;
      for (j = 0; j < dnwds; j++)
        done[j][nfree] = 0;
      cind++;
      def[nfree] = l;
      imcos[m][bcos] = nfree;
      imcos[l][nfree] = bcos;
      bcos = nfree;
      bpt[nfree] = lastd;
      fpt[lastd] = nfree;
      lastd = nfree;
      nfree = fpt[nfree];
      fpt[lastd] = 0;
    }
    else
      bcos = k;
  }
  if (fcos != bcos)
    coinc(fcos, bcos);
  WD(ccos, scanno);
  return (0);
}

int coinc(int c1, int c2)
{
  int  lc, hc, qh, qt, i, j, y, fhc, bhc, lim, him, x;
  char sw;
  lc = 1;
  while (lc != c1 && lc != c2)
    lc = fpt[lc];
  hc = lc == c1 ? c2 : c1;
  if (hc <= lo) {
    fprintf(stderr, "Impossible coincidence.\n");
    return (-1);
  }
  /* Cosets with numbers <=lo should certainly not be eliminated! */
  qh = 0;
  qt = 0;
  fhc = fpt[hc];
  bhc = bpt[hc];
  fpt[bhc] = fhc;
  if (fhc == 0)
    lastd = bhc;
  else
    bpt[fhc] = bhc;
  if (ccos == hc) {
    ccos = bhc;
    endcr = endr;
    clsd = 0;
  }
  if (lkah && lcl == hc)
    lcl = bhc;
  while (1) {
    fpt[hc] = nfree;
    nfree = hc;
    cind--;
    nelim++;
    x = 1;
    sw = 0;
    while (x <= *pno) {
      if (sw)
        i++;
      else
        i = pno[x];
      him = imcos[i][hc];
      if (him != 0) {
        j = inv[i];
        lim = imcos[i][lc];
        if (him == hc)
          him = lc;
        else
          imcos[j][him] = 0;
        if (lim == 0)
          imcos[i][lc] = him;
        else {
          if (lim == hc) {
            imcos[i][lc] = lc;
            lim = lc;
          }
          while (fpt[him] < 0)
            him = -fpt[him];
          while (fpt[lim] < 0)
            lim = -fpt[lim];
          if (him != lim) {
            y = 1;
            while (y != him && y != lim)
              y = fpt[y];
            if (y == him) {
              him = lim;
              lim = y;
            }
            if (him <= lo) {
              fprintf(stderr, "Impossible coincidence(q).\n");
              return (-1);
            }
            fhc = fpt[him];
            bhc = bpt[him];
            fpt[bhc] = fhc;
            if (fhc == 0)
              lastd = bhc;
            else
              bpt[fhc] = bhc;
            if (ccos == him) {
              ccos = bhc;
              endcr = endr;
              clsd = 0;
            }
            if (lkah && lcl == him)
              lcl = bhc;
            fpt[him] = -lim;
            if (qh == 0)
              qh = him;
            else
              bpt[qt] = him;
            qt = him;
            bpt[qt] = 0;
          }
        }
        y = imcos[i][lc];
        if (imcos[j][y] == 0)
          imcos[j][y] = lc;
      }
      if (sw)
        x++;
      sw = sw ? 0 : 1;
    }
    if (qh == 0)
      break;
    hc = qh;
    qh = bpt[qh];
    lc = -fpt[hc];
    bpt[hc] = 0;
  }
  return (0);
}
