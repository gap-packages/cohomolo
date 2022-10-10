#include "defs.h"
#include "permfns.h"

extern char cent, sym, hgst, nop, nonb[], inf1[], inf2[], inf3[], outf1[],
    outf2[];
extern short mp, mexp, mb, mnpt, prime, perm[], sv[], cp[], orb[], gbase[],
    hbase[], obase[], nbase[], lorbg[], lorbn[], lorbh[], ntno[], reg[],
    endorno[], ntorno[], tsv1[], tsv2[], tsv3[], genorb[], expcp[], fp[],
    pno[], start[], space[], ipno[], *pptr[], *svgptr[], *svhptr[], *svnptr[],
    *intorb[], *horno[], *hlorb[], *expptr[], *imorno[], *imlorb[],
    *orbperm[], *deftime[], *regsv[];
extern int sp, svsp;

extern short npt, npt1, nptd2, nbg, nbh, stg, sth, nf, nrego, nnth, lnth,
    nhorbs, ad1, hgno, lexp, *norsh, *norsim, *defim, *orlist, *imhno, *inim,
    fsth, *adpt, *bindor, *reginfo, *test;
short adno, gskno, nskno, orhct, ntct, ontct, oorhct, onf, *sva, *ap, stn,
    *tp, *itp, nnps;
char fail, bt, biggersylp, type[14];

extern FILE *ip, *op;

void err(void)
{
  fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
}

int clrop(void)
/* Clears and updates orbit permutation info */
{
  short i, j, k, l;
  for (i = 1; i <= nhorbs; i++) {
    l = *hlorb[i];
    for (j = 1; j <= l; j++) {
      k = orbperm[i][j];
      if (k != 0 && deftime[i][k] >= adno) {
        deftime[i][k] = 0;
        orbperm[i][j] = 0;
      }
    }
  }
  return (0);
}

void bind(int pno)
/* Adds a permutation in H or a newly found one in N (=N(H)) to those
   generating the orbit of gbase[gskno]. These orbit points need not be
   considered as images of new elements.
*/
{
  short i, j, k, l, m, n;
  for (i = 1; i <= npt; i++) {
    j = pptr[pno][i];
    k = bindor[i];
    l = bindor[j];
    if (k != l) {
      if (k > l) {
        m = k;
        k = l;
        l = m;
      }
      for (n = 1; n <= npt; n++)
        if (bindor[n] == l)
          bindor[n] = k;
    }
  }
}

int setskno(void)
/* Preliminary processing done for a new value of gskno. */
{
  short i, j, k, n, *p;
  /* if reg[adno]>0, we need not consider this as gskno */
  while (reg[adno] > 0) {
    adno -= 1;
    inim[gbase[adno]] = 0;
    clrop();
  }
  /* Remember initial values of orhct and oorhct. These are recalled after we
     find a new generator in N(H).
  */
  oorhct = orhct;
  ontct = ntct;
  *expcp = 1;
  gskno = adno;
  /* Set up info needed for running thro elements in search list. */
  if (sym) {
    bt = 1;
    for (i = 1; bt; i++)
      if (inim[i] == 0 && i != gbase[adno]) {
        bt = 0;
        expcp[1] = i;
      }
    bt = 1;
    for (i = npt; bt; i--)
      if (inim[i] == 0 && i != gbase[adno]) {
        bt = 0;
        start[adno + 1] = i;
      }
  }
  else {
    expcp[1] = start[adno] + 1;
    if (adno <= lexp) {
      sva = svgptr[adno];
      ap = adpt + adno;
      *ap = 1;
      while (sva[*ap] <= 0)
        (*ap)++;
      exprep(*ap, adno, sva);
    }
  }
  *pno = 0;
  if (ntorno[adno] > 0) {
    p = imlorb[orhct];
    for (j = 1; j <= *p; j++)
      p[j] = abs(p[j]);
  }
  /* Initialize the "bound" orbit, using bind, as above, and add the relevant
     elements of H to the generators of N.
   */
  for (i = 1; i <= npt; i++)
    bindor[i] = i;
  bindor[gbase[gskno]] = 0;
  n = 0;
  for (i = stn; i != -1; i = fp[i]) {
    bind(i);
    (*pno)++;
    pno[*pno] = i;
    n = i;
  }
  k = ntno[adno];
  if (k != 0) {
    for (i = sth; i != -1; i = fp[i])
      if (pptr[i][npt1] == k) {
        if (nf == -1) {
          err();
          return (-1);
        }
        for (j = 1; j <= npt + npt1; j++)
          pptr[nf][j] = pptr[i][j];
        pptr[nf][npt1] = nskno;
        (*pno)++;
        pno[*pno] = nf;
        if (stn == -1)
          stn = nf;
        else
          fp[n] = nf;
        bind(nf);
        n = nf;
        nf = fp[nf];
      }
    fp[n] = -1;
  }
  lorbn[nskno] = orbitsv(nbase[nskno], svnptr[nskno], 0);
  printf("gskno=%d, nskno=%d.\n", gskno, nskno);
  tp[npt1] = nskno;
  return (0);
}

int found(void)
/* Processing required on finding a new generator of N. */
{
  short i, j, k, *ptr, olo;
  *pno = 0;
  for (i = stn; i != -1; i = fp[i]) {
    (*pno)++;
    pno[*pno] = i;
    j = i;
  }
  if (stn == -1)
    stn = onf;
  else
    fp[j] = onf;
  (*pno)++;
  pno[*pno] = onf;
  j = (ntorno[gskno] > 0 || gskno > lnth) ? oorhct + 1 : oorhct;
  for (i = j; i <= orhct; i++) {
    ptr = imlorb[i];
    for (k = 1; k <= *ptr; k++)
      ptr[k] = abs(ptr[k]);
  }
  orhct = oorhct;
  ntct = ontct;
  fp[onf] = -1;
  bind(onf);
  if (nop) {
    tp[npt + 1] = 1;
    printvec(tp, 1);
    tp[npt + 1] = nskno;
  }
  olo = lorbn[nskno];
  lorbn[nskno] = orbitsv(nbase[nskno], svnptr[nskno], 0);
  if (prime == 0)
    biggersylp = 0;
  else
    biggersylp = ((lorbn[nskno] / olo) % prime == 0);
  onf = nf;
  if (onf == -1) {
    err();
    return (-1);
  }
  nf = fp[nf];
  tp = pptr[onf];
  itp = tp + npt1;
  for (i = 1; i < gskno; i++)
    tp[gbase[i]] = gbase[i];
  tp[npt1] = nskno;
  *expcp = 1;
  adno = gskno;
  if (sym) {
    for (i = 1; i <= npt; i++)
      inim[i] = 0;
    for (i = 1; i < adno; i++)
      inim[gbase[i]] = 1;
  }
  if (cent)
    strcpy(type, "centralizer");
  else
    strcpy(type, "normalizer");
  printf("Found element in %s. lorb=%3d.\n", type, lorbn[nskno]);
  nnps++;
  return (0);
}

void deforbperm(int bpt, int imbpt)
/* Update orbit permutation information using fact that tp[bpt]=imbpt. */
{
  short i, x, y, z;
  for (i = 1; fail == 0 && i <= orhct; i++) {
    x = horno[i][bpt];
    y = imorno[i][imbpt];
    z = orbperm[i][x];
    if (z == 0) {
      if (deftime[i][y] == 0) {
        if (hlorb[i][x] == abs(imlorb[i][y])) {
          orbperm[i][x] = y;
          deftime[i][y] = adno;
        }
        else
          fail = 1;
      }
      else
        fail = 1;
    }
    else if (z != y)
      fail = 1;
  }
}

int nprg2(void)
{
  short nbn, i, j, k, l, m, n, bpt, imbpt, ct, endo, *ptr;
  char  hornt, complete, advance, bt, orhad;
  int   quot;
  nnps = 0;
  stn = -1;
  nbn = 0;
  /* Compute a largest possible basis of N */
  for (i = 1; i <= nbg; i++)
    if (reg[i] <= 0) {
      nbn++;
      if (nbn > mb) {
        fprintf(stderr, "Too many N-base points. Increase MB.\n");
        return (-1);
      }
      nbase[nbn] = gbase[i];
    }
  printf("nbn=%d.\n", nbn);
  k = (sym) ? mb : lexp + mb;
  quot = svsp / npt - k;
  if ((sym == 0) && *obase > nbn)
    m = *obase;
  else
    m = nbn;
  if (quot < m) {
    fprintf(stderr, "svsp too small. Increase SVSP.\n");
    return (-1);
  }
  ptr = (sym) ? sv - 1 : sv + lexp * npt - 1;
  for (i = 1; i <= m + mb; i++) {
    svnptr[i] = ptr;
    ptr += npt;
  }
  if (sym)
    for (i = 1; i <= npt; i++)
      inim[i] = 0;

  /* During the main search, we attempt to construct a permutation tp, which
     we hope will lie in N(H). For each base point in turn, bpt=gbase[adno],
     we consider the image imbpt=tp[bpt], and subject it to various tests. If
     it fails any test, then the flag fail becomes true, and this image is
     known to be impossible, so we attempt to alter this image, without
     changing tp[gbase[i]] for i<adno. If impossible, we backtrack, and reduce
     adno. gskno is the G-base no for which we are currently searching in
     G[gskno], and nskno the corresponding N-base no.o orhct and ntct are
     approximately the values of orhno and ntno for the current adno. First we
     initialize everything.
  */
  j = 0;
  n = 0;
  ntct = 0;
  orhct = 0;
  for (i = 1; i < nbg; i++) {
    if (ntorno[i] > 0)
      orhct++;
    if (ntno[i] > 0)
      ntct++;
    /* regsv is used when we compute automorphism actions of tp on O, as
       described in normp1.c
    */
    if (reg[i] < 0 && ntno[i] > 0) {
      n++; /* for (j=1;j<=npt;j++) regsv[n][j]=svhptr[ntct][j]; */
      regsv[n] = svhptr[ntct];
    }
    k = gbase[i];
    for (j = 1; j <= orhct; j++) {
      l = horno[j][k];
      if (deftime[j][l] == 0) {
        orbperm[j][l] = l;
        deftime[j][l] = i;
      }
    }
    inim[k] = 1;
  }
  if (ntorno[nbg] > 0 || cent) {
    orhct++;
    ntct++;
  }
  adno = nbg;
  nskno = nbn;
  onf = nf;
  nf = fp[nf];
  if (onf == -1) {
    err();
    return (-1);
  }
  tp = pptr[onf];
  itp = tp + npt1;
  for (i = 1; i < nbg; i++)
    tp[gbase[i]] = gbase[i];
  tp[npt1] = nskno;

  printf("\nBeginning main search.\n");
  if (setskno() == -1)
    return (-1);
  if (nop) {
    op = fopen(outf2, "w");
    fprintf(op, "                    \n");
  }
  complete = 0;
  /* The main search loop follows. */
  while (complete == 0) {
    bpt = gbase[adno];
    imbpt = (sym) ? expcp[*expcp] : expimage(bpt);
    tp[bpt] = imbpt;
    hornt = (ntorno[adno] > 0);
    fail = 0;
    advance = 0;
    orhad = 0;
    if (adno == gskno) {
      if (bindor[imbpt] == 0)
        fail = 1;
    }
    /* Fails since imbpt is already in image of N */
    else if (hornt) {
      if (imlorb[orhct][imorno[orhct][imbpt]] < 0)
        fail = 1;
    }
    /* Fails since this orbit has already been ruled out as a possible image.
     */
    if (reg[adno] > 0 && imbpt != defim[adno])
      fail = 1;
    /* defim[adno] is the value of tp[bpt] which is forced by earlier values
       of tp, and automorphism information computed earlier from reginfo.
    */
    if (fail == 0)
      deforbperm(bpt, imbpt);
    /* If orbits are not permuted, deforbperm sets fail=1 */
    if (fail == 0 && ntorno[adno] < 0)
    /* Now if ntorno[adno]<0 we check that the image of O1  has kernel at
       least as big as that of the image of O.
    */
    {
      n = -ntorno[adno] - 1;
      m = (n == 0) ? 0 : ntorno[n];
      if (m > 0) {
        l = imorno[m][imbpt];
        ct = 0;
        for (i = 1; i <= npt; i++)
          if (imorno[m][i] == l) {
            ct++;
            orlist[ct] = i;
          }
        i = sth;
        while (pptr[i][npt1] > ntct) {
          for (j = 1; j <= ct; j++) {
            k = orlist[j];
            if (pptr[i][k] != k)
              fail = 1;
          }
          i = fp[i];
        }
      }
    }
    if (fail == 0 && hornt)
    /* In case adno corresponds to an H-base point, we now change the
       hbase[adno] to imbpt, and compute stabilizer orbits, as imorno,imlorb.
       If these do not correspond to the originals in H, then fail becomes
       true. Eventually they will be the images of the H orbits.
    */
    {
      if (intbase(imbpt, ntct, &sth, &nbh, hbase, lorbh, svhptr) == -1)
        return (-1);
      if (orhct < nhorbs) {
        orhct++;
        orhad = 1;
        *pno = 0;
        i = sth;
        while (pptr[i][npt1] > ntct) {
          (*pno)++;
          pno[*pno] = i;
          i = fp[i];
        }
        allorbs(test, imorno[orhct]);
        ptr = hlorb[orhct];
        if (*ptr != *test) {
          fail = 1;
          orhct--;
        }
        else {
          for (j = 1; j <= nptd2; j++) {
            norsh[j] = 0;
            norsim[j] = 0;
          }
          l = 0;
          m = 0;
          for (j = 1; j <= *ptr; j++) {
            if (*(ptr + j) > nptd2)
              l = *(ptr + j);
            else
              norsh[*(ptr + j)]++;
            if (*(test + j) > nptd2)
              m = *(test + j);
            else
              norsim[*(test + j)]++;
          }
          fail = (l != m);
          for (i = 1; fail == 0 && i <= nptd2; i++)
            fail = (norsh[i] != norsim[i]);
          if (fail)
            orhct--;
          else
            for (i = 0; i <= *test; i++)
              imlorb[orhct][i] = test[i];
        }
      }
      if (fail == 0 && reg[adno] < 0)
      /* If reg[adno]<0, we compute the corresponding F in the image, which
         must have the same size as the original F.
      */
      {
        n = ntorno[adno];
        m = imorno[n][imbpt];
        ct = 0;
        l = -reg[adno + 1];
        bt = (n == nhorbs);
        for (i = 1; i <= npt; i++)
          if (imorno[n][i] == m)
            if (bt || imlorb[n + 1][imorno[n + 1][i]] == 1) {
              ct++;
              orlist[ct] = i;
            }
        if (ct != l)
          fail = 1;
        /*   *pno=0; bt= (fail==0); i=sth;
             while (bt)
             { if (i== -1) bt=0; else
               { j=pptr[i][npt1]; if (j>=ntct) {(*pno)++; pno[*pno]=i; }
                 else if (j<ntct) bt=0; i=fp[i];
               }
             }
             if (fail==0) { i= -reg[adno]; orbitsv(imbpt,regsv[i],0); }
         */
        if (fail && orhad)
          orhct--;
      }
    }
    if (fail == 0 && (endo = endorno[adno]) > 0)
    /* if adno completes the orbit O, we compute the images of the relevant
       section of H on O, using the copies of the generators of H that we made
       at the end of nprg1(). These will be used to compute the images on
       corresponding orbits O1.
    */
    {
      if (sym == 0)
        for (i = 1; i <= npt; i++)
          tp[i] = expimage(i);
      i = fsth;
      while (pptr[i][npt1] > ntct)
        i = fp[i];
      while (i != -1) {
        *cp = 0;
        for (j = endo; fail == 0 && j <= adno; j++) {
          l = ntno[j];
          m = image(tp[pptr[i][gbase[j]]]);
          if (l > 0) {
            if (svhptr[l][m] != 0)
              addsv(m, svhptr[l]);
            else
              fail = 1;
          }
          else if (m != tp[gbase[j]])
            fail = 1;
        }
        if (fail)
          i = -1;
        else {
          k = i + 1;
          for (j = 1; j <= npt; j++)
            pptr[k][j] = backimage(j);
          i = fp[i];
          if (i != -1 && pptr[i][npt1] < ntno[endo])
            i = -1;
        }
      }
      if (fail && orhad)
        orhct--;
    }
    if (fail == 0)
    /* tp has passed all tests for this adno, so we increase adno. If
       adno=nbg, we are ready to compute tp in full, and see if it really lies
       in N(H).
    */
    {
      inim[imbpt] = 1;
      adno++;
      if (adno <= nbg) {
        advance = 1;
        if (reg[adno] > 0)
        /* if reg[adno]>0 we can compute defim[gbase[adno]] using info already
           computed, and stored in reginfo.
        */
        {
          n = reg[adno];
          m = reginfo[n];
          if (m > 0)
          /* Case 1:  gbase[adno] in set F in orbit O.  */
          {
            l = -reg[m];
            n++;
            ct = reginfo[n];
            *cp = 0;
            for (i = 1; i <= ct; i++) {
              n++;
              addsv(tp[reginfo[n]], regsv[l]);
            }
            defim[adno] = image(tp[gbase[m]]);
          }
          else
          /* Case 2:  gbase[adno] in an orbit O1.  */
          {
            n++;
            *cp = reginfo[n];
            for (k = 1; k <= *cp; k++)
              cp[k] = reginfo[n + (*cp) + 1 - k];
            defim[adno] = image(tp[-m]);
          }
        }
        if (sym) {
          bt = 1;
          for (i = npt; bt; i--)
            if (inim[i] == 0) {
              start[adno + 1] = i;
              bt = 0;
            }
          bt = 1;
          if (reg[adno] > 0) {
            k = defim[adno];
            if (inim[k]) {
              adno--;
              if (orhad)
                orhct--;
              advance = 0;
              fail = 1;
            }
            else {
              (*expcp)++;
              expcp[*expcp] = k;
            }
          }
          else
            for (i = 1; bt; i++)
              if (inim[i] == 0) {
                (*expcp)++;
                expcp[*expcp] = i;
                bt = 0;
              }
        }
      }
      else
      /* Compute tp, checking with deforbperm as we go.  */
      {
        for (i = 1; fail == 0 && i <= npt; i++)
          if (nonb[i]) {
            if (sym) {
              for (k = 1; k <= npt; k++)
                if (inim[k] == 0) {
                  j = k;
                  tp[i] = k;
                }
            }
            else {
              j = expimage(i);
              tp[i] = j;
            }
            deforbperm(i, j);
          }
        if (fail == 0) {
          invert(tp, itp);
          /* See if it lies in N(H)  (or C(H))  */
          i = (hgst) ? 0 : sth;
          while (i != -1 && fail == 0) {
            if (cent) {
              for (j = 1; fail == 0 && j <= npt; j++)
                if (tp[pptr[i][itp[j]]] != pptr[i][j])
                  fail = 1;
            }
            else {
              *cp = 0;
              for (j = 1; fail == 0 && j <= nbg; j++) {
                k = image(tp[pptr[i][gbase[j]]]);
                l = ntno[j];
                if (l != 0) {
                  if (svhptr[l][k] != 0)
                    addsv(k, svhptr[l]);
                  else
                    fail = 1;
                }
                else if (k != tp[gbase[j]])
                  fail = 1;
              }
            }
            if (hgst) {
              i++;
              if (i == hgno)
                i = -1;
            }
            else
              i = fp[i];
          }
        }
        if (fail)
          adno--;
        else {
          if (found() == -1)
            return (-1);
          if (biggersylp) {
            ad1 = gskno;
            break;
          }
        }
      }
    }
    if (advance)
    /* Change tp, changing image at current value of adno if possible. */
    {
      if (ntno[adno] > 0)
        ntct++;
      if (adno <= lexp)
        adpt[adno] = 0;
    }
    else {
      bt = 1;
      while (bt) {
        imbpt = tp[gbase[adno]];
        if (fail) {
          m = ntorno[adno];
          if (m > 0)
            n = imorno[orhct][imbpt];
        }
        if (reg[adno] <= 0 || defim[adno] != imbpt) {
          if (sym == 0) {
            if (expcp[*expcp] <= start[adno]) {
              (*expcp)++;
              expcp[*expcp] = start[adno] + 1;
              bt = 0;
            }
            else if (adno > lexp && expcp[*expcp] < start[adno + 1]) {
              expcp[*expcp]++;
              bt = 0;
            }
            if (adno <= lexp) {
              sva = svgptr[adno];
              ap = adpt + adno;
              (*ap)++;
              while (sva[*ap] <= 0 && *ap <= npt)
                (*ap)++;
              if (*ap <= npt) {
                bt = 0;
                exprep(*ap, adno, sva);
              }
              else
                bt = 1;
            }
          }
          else if (expcp[*expcp] < start[adno + 1]) {
            j = expcp[*expcp];
            inim[j] = 0;
            for (i = j + 1; bt; i++)
              if (inim[i] == 0 && (adno != gskno || i != gbase[adno])) {
                expcp[*expcp] = i;
                bt = 0;
              }
          }
        }
        if (bt) {
          j = ntorno[adno];
          if (j > 0) {
            ptr = imlorb[orhct];
            for (k = 1; k <= *ptr; k++)
              ptr[k] = abs(ptr[k]);
            ntct--;
          }
          if (sym) {
            inim[expcp[*expcp]] = 0;
            (*expcp)--;
          }
          else if (expcp[*expcp] > start[adno])
            (*expcp)--;
          adno--;
          if (adno == 0)
            bt = 0;
          else {
            if (ntorno[adno] > 0 && adno != lnth)
              orhct--;
            if (*expcp == 0)
              bt = 0;
          }
        }
      }
      clrop();
      if (*expcp == 0) {
        if (adno < ad1) {
          complete = 1;
          printf("Algorithm complete.\n");
        }
        else {
          inim[gbase[adno]] = 0;
          nskno--;
          if (setskno() == -1)
            return (-1);
        }
      }
      else if (fail) {
        if (adno == gskno) {
          k = bindor[imbpt];
          for (i = 1; i <= npt; i++)
            if (bindor[i] == k)
              bindor[i] = 0;
        }
        else if (m > 0)
          imlorb[m][n] = -abs(imlorb[m][n]);
      }
    }
  } /* Main search loop  */

  /* End of search. Add any remaining generators to N, change back to old
     base if necessary, and output.
  */
  if (cent == 0) {
    *pno = 0;
    for (i = stn; i != -1; i = fp[i]) {
      (*pno)++;
      pno[*pno] = i;
      n = i;
    }
    for (m = ad1 - 1; m >= 1; m--) {
      k = ntno[m];
      if (k != 0) {
        for (i = sth; i != -1; i = fp[i])
          if (pptr[i][npt1] == k) {
            if (nf == -1) {
              err();
              return (-1);
            }
            for (j = 1; j <= npt + npt1; j++)
              pptr[nf][j] = pptr[i][j];
            pptr[nf][npt1] = m;
            (*pno)++;
            pno[*pno] = nf;
            if (stn == -1)
              stn = nf;
            else
              fp[n] = nf;
            n = nf;
            nf = fp[nf];
          }
        fp[n] = -1;
      }
      lorbn[m] = orbitsv(nbase[m], svnptr[m], 0);
    }
  }
  if (nop) {
    fclose(op);
    if (nnps == 0)
      unlink(outf2);
    else {
      op = fopen(outf2, "r+");
      fprintf(op, "%4d%4d%4d%4d", npt, nnps, 0, 0);
      fclose(op);
    }
  }
  if (sym == 0)
    for (i = 1; i <= *obase; i++)
      if (intbase(obase[i], i, &stn, &nbn, nbase, lorbn, svnptr) == -1)
        return (-1);
  bt = 1;
  for (i = 1; i <= nbn; i++)
    if (lorbg[i] != lorbn[i]) {
      bt = 0;
      break;
    }
  *pno = 0;
  for (i = stn; i != -1; i = fp[i]) {
    (*pno)++;
    pno[*pno] = i;
  }
  if (*pno == 0) {
    printf("Group is trivial.\n");
    return (-1);
  }
  op = fopen(outf1, "w");
  fprintf(op, "%4d%4d%4d%4d\n", npt, *pno, nbn, 3);
  printbaselo(nbn, nbase, lorbn);
  printpsv(nbn, pno, svnptr);
  fclose(op);
  return (0);
}
