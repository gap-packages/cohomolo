#include "defs.h"
#include "permfns.h"

extern short npt, nf, cp[], orb[], pno[], fp[], tsv1[], tsv2[], tsv3[],
    orep[], *pptr[], mb;
short cpno;

int addperm(void)
/* This computes the perm cp as perm no cpno, and adds it to the list pno.
   Externals: nf,cpno,fp,pno,pptr.
*/
{
  short k, m, n;
  if (nf == -1) {
    fprintf(stderr, "Out of space. Increase PSP.\n");
    return (-1);
  }
  cpno = nf;
  nf = fp[nf];
  (*pno)++;
  pno[*pno] = cpno;
  m = cpno + 1;
  for (n = 1; n <= npt; n++) {
    k = image(n);
    pptr[cpno][n] = k;
    pptr[m][k] = n;
  }
  return (0);
}

int intbase(int      pt,
            int      pos,
            short *  stad,
            short *  nbad,
            short *  b,
            short *  lorb,
            short ** svptr)
/* This changes base no pos to pt (unless it is already an earlier base point.
   b is the base, lorb, svptr as usual. nbad is tha address of nb (no of
   base points). It is assumed that the perm nos are linked by the extern
   array fp, starting at perm no *stad.
   Externals:npt,mb,tsv1,tsv2,tsv3,pno,pptr,svptr,orb,orep,nf,cp.
*/
{
  char  recalc, sk1, sk2;
  short i, j, k, l, m, n, v, z, b1, b2, lo, lo1, lo2, nlo1, nlo2, np1, np2,
      fpno, lpno, nr, stpos, npt1;
  npt1 = npt + 1;
  recalc = 0;
  stpos = pos - 1;
  if (pos == *nbad + 1)
    stpos = 0;
  else
    for (i = 1; i <= *nbad; i++) {
      if (i == pos)
        stpos = 0;
      if (b[i] == pt) {
        if (stpos != 0) {
          fprintf(stderr, "Base element there.\n");
          return (1);
        }
        stpos = i;
      }
    }

  /* If stpos=0, then bno is not at present a base point, so we introduce it
     as base[nb+1], and put stpos=nb+1. Otherwise, it is already present as
     base[stpos]. The idea is to swap base points l,l+1  for l=stpos-1,...pos
     successively.
  */
  if (stpos == 0) {
    (*nbad)++;
    if (*nbad > mb) {
      fprintf(stderr, "nb too big. Increase SVSP (or MB).\n");
      return (-1);
    }
    stpos = *nbad;
    b[*nbad] = pt;
    lorb[*nbad] = 1;
    for (n = 1; n <= npt; n++)
      svptr[*nbad][n] = 0;
    svptr[*nbad][pt] = -1;
  }
  for (l = stpos - 1; l >= pos; l--) {
    b1 = b[l];
    b2 = b[l + 1];
    lo1 = lorb[l];
    lo2 = lorb[l + 1];
    b[l] = b2;
    b[l + 1] = b1;
    /* Now we list the permutations that may need to be changed */
    *pno = 0;
    sk1 = 1;
    sk2 = 1;
    for (i = *stad; i != -1; i = fp[i]) {
      k = pptr[i][npt1];
      if (sk1 && k <= l + 1) {
        np1 = *pno;
        sk1 = 0;
      }
      if (sk2 && k <= l) {
        np2 = *pno;
        sk2 = 0;
      }
      if (k < l)
        break;
      (*pno)++;
      pno[*pno] = i;
    }
    if (sk1)
      np1 = *pno;
    if (sk2)
      np2 = *pno;
    if (lo1 == 1) {
      for (n = 1; n <= npt; n++) {
        svptr[l][n] = svptr[l + 1][n];
        svptr[l + 1][n] = 0;
      }
      svptr[l + 1][b1] = -1;
      lorb[l] = lo2;
      lorb[l + 1] = 1;
      for (i = np1 + 1; i <= np2; i++)
        pptr[pno[i]][npt1] = l;
      continue;
    }
    lo = orbitsv(b2, tsv1, 0);
    if (lo == 1) {
      for (n = 1; n <= npt; n++) {
        svptr[l + 1][n] = svptr[l][n];
        svptr[l][n] = 0;
      }
      svptr[l][b2] = -1;
      lorb[l] = 1;
      lorb[l + 1] = lo1;
      for (i = np2 + 1; i <= *pno; i++)
        pptr[pno[i]][npt1] = l + 1;
      continue;
    }
    nlo1 = lo;
    lorb[l] = lo;
    nlo2 = lo1 * lo2 / lo;
    lorb[l + 1] = nlo2;
    if (nlo2 == 1) {
      for (n = 1; n <= npt; n++) {
        svptr[l][n] = tsv1[n];
        svptr[l + 1][n] = 0;
      }
      svptr[l + 1][b1] = -1;
      for (i = np1 + 1; i <= np2; i++)
        pptr[pno[i]][npt1] = l;
      continue;
    }
    /* The easy cases have now been handled! */
    if (recalc)
      orbitsv(b1, svptr[l], 0);
    else
      recalc = 1;
    lpno = pno[*pno];
    *pno = np2;
    orbitsv(b2, tsv2, 0);
    *pno = np1;
    for (n = 1; n <= npt; n++)
      tsv3[n] = 1;
    nr = 0;
    for (n = 1; n <= npt; n++)
      if (svptr[l][n] > 0 && tsv3[n]) {
        nr++;
        orb[1] = n;
        orep[nr] = n;
        tsv3[n] = 0;
        m = 1;
        for (j = 1; j <= m; j++) {
          z = orb[j];
          for (k = 1; k <= *pno; k++) {
            v = pptr[pno[k]][z];
            if (tsv3[v]) {
              tsv3[v] = 0;
              m++;
              orb[m] = v;
            }
          }
        }
      }
    for (n = 1; n <= npt; n++)
      tsv3[n] = 0;
    tsv3[b1] = -1;
    lo = 1;
    orb[1] = b1;
    fpno = -1;
    for (i = 1; i <= nr; i++) {
      *cp = 0;
      addsv(orep[i], svptr[l]);
      j = image(b2);
      if (tsv2[j] != 0) {
        addsv(j, tsv2);
        j = image(b1);
        if (tsv3[j] == 0) {
          if (fpno == -1)
            fpno = nf;
          else
            fp[cpno] = nf;
          if ((addperm()) == -1)
            return (-1);
          pptr[cpno][npt1] = l + 1;
          if ((lo = orbitsv(b1, tsv3, lo)) == nlo2)
            break;
        }
      }
    }
    for (n = 1; n <= npt; n++) {
      svptr[l + 1][n] = tsv3[n];
      tsv3[n] = 0;
    }
    tsv3[b2] = -1;
    lo = 1;
    orb[1] = b2;
    for (n = 1; n <= npt; n++)
      if (tsv1[n] > 0 && tsv3[n] == 0) {
        *cp = 0;
        addsv(n, tsv1);
        fp[cpno] = nf;
        if (addperm() == -1)
          return (-1);
        pptr[cpno][npt1] = l;
        if ((lo = orbitsv(b2, tsv3, lo)) == nlo1)
          break;
      }
    for (n = 1; n <= npt; n++)
      svptr[l][n] = tsv3[n];
    fp[cpno] = fp[lpno];
    fp[lpno] = nf;
    if (np1 == 0) {
      nf = *stad;
      *stad = fpno;
    }
    else {
      j = pno[np1];
      nf = fp[j];
      fp[j] = fpno;
    }
  }
  while (lorb[*nbad] == 1 && b[*nbad] != pt)
    (*nbad)--;
  if (recalc && pos > 1) {
    *pno = 0;
    j = pos - 1;
    for (i = *stad; i != -1; i = fp[i]) {
      k = pptr[i][npt1];
      if (k < j) {
        orbitsv(b[j], svptr[j], 0);
        j = k;
      }
      (*pno)++;
      pno[*pno] = i;
    }
    orbitsv(b[j], svptr[j], 0);
  }
  return (0);
}
