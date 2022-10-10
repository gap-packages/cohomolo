#include "defs.h"

#define tmalloc(D, T, N)                                                     \
  {                                                                          \
    D = (T *)malloc(sizeof(T) * (N));                                        \
    if (D == 0) {                                                            \
      fprintf(stderr, "Out of space.\n");                                    \
      return (-1);                                                           \
    }                                                                        \
  }
#define tfree(D)                                                             \
  {                                                                          \
    if (D)                                                                   \
      free((char *)D);                                                       \
    D = 0;                                                                   \
  }

extern char  mult;
extern short wd1[], wd2[], wd3[], wd4[], *cst, *cend, ***coeff, fullsc, clsd,
    lkah, conch, **cco, **def;
extern short invg[], mwdl, endr, *fpt, *bpt, nelim, lo, **imcos, dim, ccos,
    lastd, cind, nfree, stcr, endcr, fcos, bcos, lcl, np2, *rel, ng;

int scanrel(void)
{
  short i, j, k, l, m, necr, *s, *t;
  short compfsc, *p, *q, esc;
  fullsc = 1;
  stcr = endcr + 2;
  endcr += (1 + rel[stcr - 1]);
  esc = 0;
  *wd1 = 0;
  *wd2 = 0;
  p = wd1 + mwdl;
  q = p + dim;
  while (++p <= q)
    *p = 0;
  p = wd2 + mwdl;
  q = p + dim;
  while (++p <= q)
    *p = 0;
  p = wd2 + mwdl;
  s = rel + endcr + 1;
  t = s + *s;
  necr = endcr + *s + 1;
  while (++s < t) {
    p[*s] = *(s + 1);
    s++;
  }
  if (mult == 0 && ccos > 1)
    action(p, def[ccos]);
  fcos = ccos;
  bcos = ccos;
  compfsc = 1;
  for (i = stcr; i <= endcr; i++) {
    k = imcos[rel[i]][fcos];
    if (k == 0) {
      compfsc = 0;
      break;
    }
    if (concat(wd1, coeff[rel[i]][fcos]) == -1)
      return (-1);
    fcos = k;
  }
  if (compfsc) {
    if (fcos != bcos)
      esc = coinc(fcos, bcos);
    else if (conch)
      if (concheck() == -1)
        return (-1);
    if (esc == -1)
      return (-1);
    if (esc == 0)
      endcr = necr;
    return (0);
  }
  stcr = i;
  for (i = endcr; i >= stcr; i--) {
    l = rel[i];
    m = invg[l];
    k = imcos[m][bcos];
    if (k == 0) {
      if (i == stcr) {
        k = imcos[l][fcos];
        if (k == 0) {
          imcos[l][fcos] = bcos;
          imcos[m][bcos] = fcos;
          if (ainvb(wd1, wd2, wd3) == -1)
            return (-1);
          if (*wd3 > 2)
            reduce(wd3);
          if (comp(wd3, coeff[l] + fcos) == -1)
            return (-1);
          invwd(wd3, wd2);
          if (comp(wd2, coeff[m] + bcos) == -1)
            return (-1);
        }
        else {
          if (concat(wd1, coeff[l][fcos]) == -1)
            return (-1);
          if (k != bcos)
            esc = coinc(k, bcos);
          else if (conch)
            if (concheck() == -1)
              return (-1);
          if (esc == -1)
            return (-1);
        }
        if (esc == 0)
          endcr = necr;
        return (0);
      }
      if (lkah || nfree == 0) {
        fullsc = 0;
        if (esc == 0)
          endcr = necr;
        return (0);
      }
      for (j = 0; j < ng; j++)
        imcos[j][nfree] = 0;
      cind++;
      if (mult == 0) {
        tmalloc(def[nfree], short, def[bcos][0] + 2);
        def[nfree][0] = def[bcos][0] + 1;
        def[nfree][1] = l;
        for (j = 1; j <= def[bcos][0]; j++)
          def[nfree][j + 1] = def[bcos][j];
      }
      imcos[m][bcos] = nfree;
      imcos[l][nfree] = bcos;
      bcos = nfree;
      bpt[nfree] = lastd;
      fpt[lastd] = nfree;
      lastd = nfree;
      nfree = fpt[nfree];
      fpt[lastd] = 0;
    }
    else {
      if (concat(wd2, coeff[m][bcos]) == -1)
        return (-1);
      bcos = k;
    }
  }
  if (fcos != bcos)
    esc = coinc(fcos, bcos);
  else if (conch)
    if (concheck() == -1)
      return (-1);
  if (esc == -1)
    return (-1);
  if (esc == 0)
    endcr = necr;
  return (0);
}

int coinc(int c1, int c2)
{
  short  lc, hc, qh, qt, i, j, y, fhc, bhc, lim, him, ret;
  short *ocend, esc;
  esc = 0;
  ocend = cend;
  lc = 1;
  while (lc != c1 && lc != c2)
    lc = fpt[lc];
  if (lc == c1) {
    hc = c2;
    ret = ainvb(wd2, wd1, wd3);
  }
  else {
    hc = c1;
    ret = ainvb(wd1, wd2, wd3);
  }
  if (ret == -1)
    return (-1);
  if (hc <= lo) {
    fprintf(stderr, "Impossible coincidence.\n");
    return (-1);
  }
  if (compb(wd3, cco + hc) == -1)
    return (-1);
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
    esc = 1;
    clsd = 0;
  }
  if (lkah && lcl == hc)
    lcl = bhc;
  while (1) {
    fpt[hc] = nfree;
    nfree = hc;
    cind--;
    nelim++;
    for (i = 0; i < ng; i++) {
      him = imcos[i][hc];
      if (him != 0) {
        expand(wd1, coeff[i][hc]);
        j = invg[i];
        lim = imcos[i][lc];
        if (him == hc) {
          him = lc;
          if (concat(wd1, cco[hc]) == -1)
            return (-1);
        }
        else
          imcos[j][him] = 0;
        if (lim == 0) {
          imcos[i][lc] = him;
          if (ainvb(wd3, wd1, wd2) == -1)
            return (-1);
          if (*wd2 > 2)
            reduce(wd2);
          if (comp(wd2, coeff[i] + lc) == -1)
            return (-1);
        }
        else {
          expand(wd2, coeff[i][lc]);
          if (lim == hc) {
            imcos[i][lc] = lc;
            lim = lc;
            if (concat(wd2, cco[hc]) == -1)
              return (-1);
            if (*wd2 > 2)
              reduce(wd2);
            if (comp(wd2, coeff[i] + lc) == -1)
              return (-1);
          }
          while (fpt[him] < 0) {
            if (concat(wd1, cco[him]) == -1)
              return (-1);
            him = -fpt[him];
          }
          while (fpt[lim] < 0) {
            if (concat(wd2, cco[lim]) == -1)
              return (-1);
            lim = -fpt[lim];
          }
          if (him != lim) {
            if (ainvb(wd3, wd1, wd4) == -1)
              return (-1);
            y = 1;
            while (y != him && y != lim)
              y = fpt[y];
            if (y == him) {
              him = lim;
              lim = y;
              ret = ainvb(wd2, wd4, wd1);
            }
            else
              ret = ainvb(wd4, wd2, wd1);
            if (ret == -1)
              return (-1);
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
              esc = 1;
              clsd = 0;
            }
            if (lkah && lcl == him)
              lcl = bhc;
            if (*wd1 > 2)
              reduce(wd1);
            if (compb(wd1, cco + him) == -1)
              return (-1);
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
        if (imcos[j][y] == 0) {
          imcos[j][y] = lc;
          expand(wd1, coeff[i][lc]);
          invwd(wd1, wd2);
          if (comp(wd2, coeff[j] + y) == -1)
            return (-1);
        }
      }
      coeff[i][hc] = 0;
    }
    if (qh == 0)
      break;
    hc = qh;
    qh = bpt[qh];
    lc = -fpt[hc];
    bpt[hc] = 0;
    expand(wd3, cco[hc]);
  }
  cend = ocend;
  return (esc);
}
