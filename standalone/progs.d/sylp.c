#include "defs.h"
#include "permfns.h"

extern char  inf[], outf[];
extern short mp, mxexp, mb, mnpt, perm[], svp2[], svg2[], cp[], orb[],
    expperm[], adpt[], expcp[], lorbg[], lorbp[], lorbdef[], base[], ntfpt[],
    ntbpt[], start[], invbase[], pno[], facord[], orno[], lporb[], deftime[],
    orbperm[], genorb[], *pptr[], *svgptr[], *svpptr[], *expptr[];
extern int psp, expsp, svsp;

short npt, p, im, cb, adno, *tp;
char  ok;
FILE *ip, *op;
/* This program and normrun involve backtrack searches. These can be speeded
   up by storing the perms for each coset rep in the stabilizer chain, which
   can be computed using Schreier vectors. As many such perms as space allows
   are stored in this way, in array expperm, using perm pointer expptr, and
   expcp corresponds to cp.
*/

int sylprog(void)
{
  char  nontriv, bt, seek, b, incadno;
  short i, j, k, l, lnt, fnt, mxp, mexp, nperms, nb, ct, np2, bno, stp, lexp,
      *z, *itp, *ap, *sva;
  int quot;
  if ((ip = fopen(inf, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf);
    return (-1);
  }
  printf("Input prime!    ");
  scanf("%hd", &p);
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &l);
  if (npt > mnpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    return (-1);
  }
  if (nb > mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (nb * npt > svsp) {
    fprintf(stderr, "Not enough sv space.Increase SVSP.\n");
    return (-1);
  }
  if (l <= 0) {
    fprintf(stderr, "Wrong input format.\n");
    return (-1);
  }
  readbaselo(nb, base, lorbg);
  for (i = 1; i <= npt; i++)
    invbase[i] = 0;
  for (i = 1; i <= nb; i++)
    invbase[base[i]] = i;
  lnt = 0;
  fnt = 0;
  /* Determine order of sylp subgroup . Links ntfpt and ntbpt are used to
     bypass trivial indices in the stab chain.
  */
  for (i = nb; i >= 1; i--) {
    j = lorbg[i];
    lorbp[i] = 1;
    if (j > 1) {
      if (lnt == 0) {
        lnt = i;
        k = i;
        ntfpt[lnt] = nb + 1;
      }
      else {
        ntbpt[k] = i;
        ntfpt[i] = k;
        k = i;
      }
      while (j % p == 0) {
        lorbp[i] *= p;
        j /= p;
      }
      if (lorbp[i] > 1)
        fnt = i;
    }
    else
      invbase[base[i]] = 0;
  }
  ntbpt[fnt] = 0;
  if (fnt == 0) {
    fprintf(stderr, "Sylp-group is trivial!\n");
    return (-1);
  }
  printf("Orbit lengths of Sylp-group:\n");
  for (i = 1; i <= nb; i++)
    printf("%4d", lorbp[i]);
  printf(".\n");

  quot = psp / (npt + 1);
  if (quot > mp)
    quot = mp;
  mxp = quot;
  quot = expsp / npt;
  if (quot > mxexp)
    quot = mxexp;
  mexp = quot;
  np2 = 2 * nperms - 1;
  if (np2 >= mxp) {
    fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
    return (-1);
  }
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + (npt + 1) * i - 1;
  for (i = 1; i <= nb; i++)
    svgptr[i] = svg2 + npt * (i - 1) - 1;
  for (i = 1; i <= nb; i++)
    svpptr[i] = svp2 + npt * (i - 1) - 1;
  readpsv(0, nb, nperms, svgptr);
  fclose(ip);

  /* We now determine how much space we have to store cosetrep perms, and
     calculate these with the function exprep
  */
  lexp = 0;
  ct = 0;
  for (i = nb; i > fnt; i--) {
    ct += (lorbg[i] - 1);
    if (ct + i >= mexp) {
      lexp = i;
      break;
    }
  }
  if (lexp == nb) {
    fprintf(stderr, "Not enough expanding space. Increase MXEXP.\n");
    return (-1);
  }
  if (lexp == 0)
    lexp = fnt;
  z = expperm - 1;
  i = 0;
  while (i != lexp) {
    i = (i == 0) ? fnt : ntfpt[i];
    expptr[i] = z;
    z += npt;
    start[i] = i - 1;
  }
  printf("lexp=%d.\n", lexp);
  start[lexp + 1] = lexp;
  k = lexp;
  for (i = lexp + 1; i <= lnt; i++) {
    l = lorbg[i];
    start[i + 1] = start[i] + l - 1;
    if (l > 1)
      for (j = 1; j <= npt; j++)
        if (svgptr[i][j] > 0) {
          k++;
          expptr[k] = z;
          z += npt;
          exprep(j, k, svgptr[i]);
        }
  }
  bt = 0;
  nontriv = 0;
  *pno = 0;
  deftime[0] = 0;
  stp = np2 + 1;
  for (i = 1; i <= nb; i++) {
    for (j = 1; j <= npt; j++) {
      svpptr[i][j] = 0;
      svpptr[i][base[i]] = -1;
    }
    lorbdef[i] = 1;
  }

  /* Starting search. lorbdef records the length of the orbit ofthe Sylp group
     found so far. tp is a perm that we are currently considering as a
     possible new element of Sylp. nontriv is set true as soon as the first
     element of Sylp has been found. adno is then number in the stabilizer
     chain for which the coset reps are currently being advanced in the
     search.
  */
  for (bno = lnt; bno >= fnt; bno = ntbpt[bno]) {
    printf("bno=%d.\n", bno);
    orb[1] = base[bno];
    adno = bno;
    while (lorbdef[bno] < lorbp[bno]) {
      *expcp = 0;
      tp = pptr[np2 + 1];
      itp = pptr[np2 + 2];
      if (np2 + 2 >= mxp) {
        fprintf(stderr, "Out of space. Increase PSP.\n");
        return (-1);
      }
      for (i = 1; i <= npt; i++)
        tp[i] = i;
      tp[npt + 1] = bno;

      /* If nontriv, then we compute all orbits of Sylp. The new element must
         permute these orbits.
      */
      if (nontriv) {
        allorbs(lporb, orno);
        bt = 0;
        for (i = 1; i <= *lporb; i++) {
          orbperm[i] = 0;
          deftime[i] = 0;
        }
      }
      seek = 1;
      while (seek) {
        for (i = 1; i <= npt; i++)
          itp[i] = 0;
        b = 1;
        ok = 1;
        while (1) {
          k = expcp[*expcp];
          /* Advance to next element in search. ok is set false if this
             element is eliminated from membership of Sylp, and we advance
             again.
          */
          if (b && (*expcp == 0 || k <= start[adno])) {
            (*expcp)++;
            if (adno <= lexp) {
              adpt[adno] = 0;
              expcp[*expcp] = adno;
            }
            else
              expcp[*expcp] = start[adno];
            b = 0;
          }
          else {
            if (adno <= lexp) {
              sva = svgptr[adno];
              ap = adpt + adno;
              (*ap)++;
              while (sva[*ap] <= 0 && *ap <= npt)
                (*ap)++;
              if (*ap <= npt) {
                exprep(*ap, adno, sva);
                break;
              }
            }
            else if (k < start[adno + 1]) {
              expcp[*expcp]++;
              break;
            }
            if (*expcp == 1) {
              fprintf(stderr,
                      "Premature end of search. Probably %d is not prime.\n",
                      p);
              return (-1);
            }
            (*expcp)--;
            bt = 1;
            adno = ntbpt[adno];
            b = 1;
            if (adno == bno && nontriv) {
              im = expimage(base[adno]);
              l = orno[im];
              if (lporb[l] > 1)
                for (i = 1; i <= npt; i++)
                  if (orno[i] == l)
                    svpptr[adno][i] = -2;
              /* Setting sv= -2 when bno=adno avoids testing more elements
                 than necessary which map the main orbit of base[bno] to
                 another.
              */
            }
          }
        }
        incadno = 1;
        while (incadno) {
          cb = base[adno];
          im = expimage(cb);
          if (adno == bno && svpptr[adno][im] != 0) {
            ok = 0;
            break;
          }
          else {
            j = im;
            k = invbase[j];
            ct = 1;
            while (k > 0 && k < adno) {
              j = tp[j];
              k = invbase[j];
              ct++;
            }
            if (k == adno && ppower(ct) == 0) {
              ok = 0;
              break;
            }
            /* Because this element cannot be a p-element */
          }
          if (nontriv) {
            if (bt)
              for (i = 1; i <= *lporb; i++) {
                j = orbperm[i];
                if (deftime[j] >= adno) {
                  orbperm[i] = 0;
                  deftime[j] = 0;
                }
              }
            bt = 0;
            deforbperm();
            if (ok == 0)
              break;
            /* This means that the element does not permute the orbits */
          }
          else
            tp[cb] = im;
          if (adno == lnt) {
            incadno = 0;
            for (cb = 1; cb <= npt; cb++)
              if (invbase[cb] == 0) {
                im = expimage(cb);
                if (nontriv) {
                  deforbperm();
                  if (ok == 0) {
                    bt = 1;
                    break;
                  }
                }
                else
                  tp[cb] = im;
              }
          }
          else
            adno = ntfpt[adno];
        }
        if (ok)
          for (i = 1; i <= npt; i++)
            if (itp[i] == 0) {
              j = i;
              l = 1;
              k = tp[j];
              itp[k] = j;
              while (k != i) {
                l++;
                j = k;
                k = tp[j];
                itp[k] = j;
              }
              if (ppower(l) == 0) {
                ok = 0;
                bt = 1;
                break;
              }
            }
        /* if ok then tp is a p-element permuting the orbits. The final test
           is to check that it normalizes Sylp.
        */
        if (ok && nontriv)
          for (i = stp; i < np2 && ok; i += 2) {
            *cp = 0;
            for (j = fnt; j <= lnt; j = ntfpt[j]) {
              k = image(tp[pptr[i][itp[base[j]]]]);
              if (svpptr[j][k] != 0)
                addsv(k, svpptr[j]);
              else {
                ok = 0;
                bt = 1;
                break;
              }
            }
          }
        /* if ok then tp is the new element of Sylp */
        if (ok) {
          (*pno)++;
          pno[*pno] = np2 + 1;
          k = lorbdef[bno];
          l = orbitsv(base[bno], svpptr[bno], k);
          facord[np2 + 1] = l / k;
          lorbdef[bno] = l;
          np2 += 2;
          adno = bno;
          nontriv = 1;
          seek = 0;
          printf("Found element in Sylp-group. Lorbdef=%d.\n", l);
          for (i = 1; i <= npt; i++)
            if (svpptr[bno][i] == -2)
              svpptr[bno][i] = 0;
        }
      }
    }
    for (i = 1; i <= npt; i++)
      if (svpptr[bno][i] == -2)
        svpptr[bno][i] = 0;
  }
  op = fopen(outf, "w");
  fprintf(op, "%4d%4d%4d%4d\n", npt, *pno, nb, 4);
  printbaselo(nb, base, lorbp);
  printpsv(nb, pno, svpptr);
  fprintf(op, "%3d", p);
  for (i = stp; i < np2; i += 2)
    fprintf(op, "%4d", facord[i]);
  fprintf(op, "\n");
  return (0);
}

int ppower(int x)
/* This tests whether x is a power of p */
{
  while (x != 1) {
    if (x % p == 0)
      x /= p;
    else
      return (0);
  }
  return (1);
}

int deforbperm(void)
/* This tests whether the fact that tp sends cb (=base[adno]) to im
   contradicts permutation of the orbits. If so then ok is set false.
*/
{
  short i, j, k;
  i = orno[cb];
  j = orno[im];
  k = orbperm[i];
  if (k == 0) {
    if (deftime[j] == 0)
      if (lporb[i] == lporb[j]) {
        orbperm[i] = j;
        deftime[j] = adno;
      }
      else
        ok = 0;
    else
      ok = 0;
  }
  else if (k != j)
    ok = 0;
  if (ok)
    tp[cb] = im;
  return (0);
}
