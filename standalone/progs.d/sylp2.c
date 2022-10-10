#include "defs.h"
#include "permfns.h"

extern char  inf1[], inf2[], outf1[], chpar;
extern short mp, mexp, mb, mnpt, par1, par2, par3, par4, prime, perm[], sv[],
    cp[], orb[], space[], expcp[], lorbg[], lorbh[], gbase[], start[], pno[],
    genorb[], *pptr[], *svgptr[], *svhptr[], *expptr[], obase[], nbase[],
    ntorno[], lorbn[], reg[], ipno[], tsv1[], tsv2[], tsv3[], orep[],
    *intorb[];
extern int psp, sp, svsp;

short npt, im, cb, adno, *tp, *itp, nb, bno, lexp, np2, *adpt, *lorbdef,
    *ntfpt, *ntbpt, *invbase, *facord, *orno, *lporb, *deft, *orbp, **gorb,
    bnoorno, opno, opct;
char  ok, bt;
FILE *ip, *op;
/* This program and normrun involve backtrack searches. These can be speeded
   up by storing the perms for each coset rep in the stabilizer chain, which
   can be computed using Schreier vectors. As many such perms as space allows
   are stored in this way, in array space, using perm pointer expptr, and
   expcp corresponds to cp.
*/

int sylprog(int x)
/* x!=0 means that a subgroup of the sylp-group has already been computed, and
   is stored in inf2. The search should now start at bno=abs(x). sylnorm
   always returns the current bno, if it exits in order to compute the
   normalizer. It returns 0 if it completes the computation of the sylp-group
   If x>0, its normalizer has also been computed, and lies in outf1,
   and so the next element of the sylp-group will be sought from this
   normalizer. Any p-element will do in this case.
*/
{
  char  nontriv, seek, b, incadno, comm, pow, id;
  short i, j, k, l, m, n, lnt, fnt, mxp, mnb, mxexp, nperms, ct, stp, *z, *ap,
      *sva, ct1, ct2, ct3, ord, kord, pt, stbno, skct, *tpk, *itpk, *commp,
      *icommp, *spptr;
  int quot;
  adpt = obase;
  lorbdef = nbase;
  ntfpt = ntorno;
  ntbpt = lorbn;
  invbase = reg;
  facord = ipno;
  orno = tsv1;
  lporb = tsv2;
  deft = tsv3;
  orbp = orep, gorb = intorb;
  if (chpar)
  /* par1,par2 and par3 are used when seeking random elts. See line 152
     (approx). When searching for an element normalizing a subgroup, we give
     up after par4 attempts, and use normalizer program to compute the
     complete normalizer first.
  */
  {
    printf(
        "Choose values of par1,par2,par3,par4. (Defaults are %d %d %d %d.)\n",
        par1, par2, par3, par4);
    scanf("%hd%hd%hd%hd", &par1, &par2, &par3, &par4);
  }
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    return (-1);
  }
  if (x == 0) {
    printf("Input prime!    ");
    scanf("%hd", &prime);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &l);
  if (npt > mnpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    return (-1);
  }
  if (nb > mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (nb * npt * 2 > svsp) {
    fprintf(stderr, "Out of sv space. Increase SVSP.\n");
    return (-1);
  }
  if (l <= 2) {
    fprintf(stderr, "Wrong input format.\n");
    return (-1);
  }
  readbaselo(nb, gbase, lorbg);
  for (i = 1; i <= npt; i++)
    invbase[i] = 0;
  for (i = 1; i <= nb; i++)
    invbase[gbase[i]] = i;
  lnt = 0;
  fnt = 0;
  /* Determine order of sylp subgroup . Links ntfpt and ntbpt are used to
     bypass trivial indices in the stab chain.
  */
  for (i = nb; i >= 1; i--) {
    j = lorbg[i];
    lorbh[i] = 1;
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
      while (j % prime == 0) {
        lorbh[i] *= prime;
        j /= prime;
      }
      if (lorbh[i] > 1)
        fnt = i;
    }
    else
      invbase[gbase[i]] = 0;
  }
  ntbpt[fnt] = 0;
  ntfpt[0] = fnt;
  if (fnt == 0) {
    fprintf(stderr, "Sylp-group is trivial!\n");
    return (-1);
  }
  printf("Orbit lengths of Sylp-group:\n");
  for (i = 1; i <= nb; i++)
    printf("%4d", lorbh[i]);
  printf(".\n");

  quot = psp / (npt + 1);
  if (quot > mp)
    quot = mp;
  mxp = quot;
  quot = sp / npt;
  if (quot > mexp)
    quot = mexp;
  mxexp = quot;
  np2 = 2 * nperms - 1;
  if (np2 >= mxp) {
    fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
    return (-1);
  }
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + (npt + 1) * i - 1;
  for (i = 1; i <= nb; i++)
    svgptr[i] = sv + npt * (i - 1) - 1;
  for (i = 1; i <= nb; i++)
    svhptr[i] = sv + npt * (i + nb - 1) - 1;
  readpsv(0, nb, nperms, svgptr);
  fclose(ip);
  tpk = space - 1;
  itpk = tpk + npt;
  commp = itpk + npt;
  icommp = commp + npt;
  spptr = icommp + npt + 1;
  sp -= 4 * npt;
  for (i = fnt; i <= nb; i = ntfpt[i]) {
    gorb[i] = spptr;
    sva = svgptr[i];
    for (j = 1; j <= npt; j++)
      if (sva[j] != 0)
        *(spptr++) = j;
    sp -= lorbg[i];
  }

  /* We now determine how much space we have to store cosetrep perms, and
     calculate these with the function exprep
  */
  lexp = 0;
  ct = 0;
  for (i = nb; i > fnt; i--) {
    ct += (lorbg[i] - 1);
    if (ct + i >= mxexp) {
      lexp = i;
      break;
    }
  }
  if (lexp == 0)
    lexp = fnt;
  z = spptr - 1;
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
  /* nontriv nonzero means in this case that any p-element in G will do; i.e.
     we are not looking for an element in the normalizer of the group found so
     far. This is the case either when x=0, right at the beginning, or when
     x>0, and we have already computed this normalizer.
  */
  nontriv = (x >= 0) ? 0 : 1;
  *pno = 0;
  deft[0] = 0;
  stp = np2 + 1;
  if (x != 0) {
    ip = fopen(inf2, "r");
    fscanf(ip, "%hd%hd%hd%hd", &i, &nperms, &i, &i);
    readbaselo(nb, gbase, lorbdef);
    readpsv(stp, nb, nperms, svhptr);
    for (i = 1; i <= nperms; i++) {
      (*pno)++;
      pno[*pno] = ++np2;
      fscanf(ip, "%hd", facord + np2);
      np2++;
    }
    fclose(ip);
    k = abs(x);
    l = 0;
    for (i = 1; i <= npt; i++)
      if (svhptr[k][i] != 0)
        orb[++l] = i;
  }
  else
    for (i = 1; i <= nb; i++) {
      for (j = 1; j <= npt; j++) {
        svhptr[i][j] = 0;
        svhptr[i][gbase[i]] = -1;
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
  stbno = (x == 0) ? lnt : abs(x);
  for (bno = stbno; bno >= fnt; bno = ntbpt[bno]) {
    printf("bno=%d.\n", bno);
    if (bno != abs(x))
      orb[1] = gbase[bno];
    adno = bno;
    while (lorbdef[bno] < lorbh[bno]) {
      if (np2 + 2 >= mxp) {
        fprintf(stderr, "Out of space. Increase PSP.\n");
        return (-1);
      }
      *expcp = 0;
      tp = pptr[np2 + 1];
      itp = pptr[np2 + 2];
      tp[npt + 1] = bno;
      /* When any p-element suffices, we seek random elements until one of
         them has order a multiple of p. If after par1 elements, we have not
         found one, we try commutators of the last element with new elements.
         If after par2 such attempts we still have not found one, we try
         commutators with powers of the last element. Then if after par3 tries
         we are still unsuccessful, we try a new random element, and use
         commutators with that and so on...
      */
      if (nontriv == 0) {
        ct1 = 0;
        seek = 1;
        comm = 0;
        while (seek) {
          (*expcp) = 0;
          ranelt();
          ord = findord(tp, itp);
          if (ord % prime != 0) {
            ct1++;
            if (comm) {
              id = 1;
              for (i = 1; i <= npt; i++) {
                j = itpk[itp[tpk[tp[i]]]];
                if (id && j != i)
                  id = 0;
                commp[i] = j;
                icommp[i] = 0;
              }
              if (id) {
                ct2++;
                if (ct2 > par2)
                  comm = 0;
                continue;
              }
              ct2 = 0;
              ord = findord(commp, icommp);
              if (ord % prime != 0) {
                ct3++;
                if (ct3 > par3) {
                  if (pow) {
                    comm = 0;
                    continue;
                  }
                  printf("Trying powers.\n");
                  pow = 1;
                  ct3 = 0;
                  ct2 = 0;
                  for (i = 2; i <= kord; i++)
                    if (kord % i == 0) {
                      if (kord == i) {
                        comm = 0;
                        continue;
                      }
                      k = kord / i;
                      for (i = 1; i <= npt; i++) {
                        pt = i;
                        for (j = 1; j <= k; j++)
                          pt = tpk[pt];
                        itpk[pt] = i;
                      }
                      invert(itpk, tpk);
                    }
                }
              }
              else {
                for (i = 1; i <= npt; i++) {
                  tp[i] = commp[i];
                  itp[i] = icommp[i];
                }
                seek = 0;
              }
            } /* if (comm) */
            else if (ct1 > par1 && ord > 1) {
              printf("Trying commutators.\n");
              ct2 = 0;
              ct3 = 0;
              kord = ord;
              comm = 1;
              pow = 0;
              for (i = 1; i <= npt; i++) {
                tpk[i] = tp[i];
                itpk[i] = itp[i];
              }
            }
          } /* if (ord%prime... */
          else
            seek = 0;
          if (seek == 0) {
            while (ord % prime == 0)
              ord /= prime;
            if (ord > 1) {
              for (i = 1; i <= npt; i++) {
                pt = i;
                for (j = 1; j <= ord; j++)
                  pt = tp[pt];
                itp[pt] = i;
              }
              invert(itp, tp);
            }
            if (svhptr[bno][tp[gbase[bno]]] != 0) {
              seek = 1;
              comm = 0;
            }
            else {
              nontriv = 1;
              fndelt();
            }
          }
        } /* while (seek) */
        continue;
      } /* if (nontriv==0) */
      for (i = 1; i <= npt; i++)
        tp[i] = i;

      /* Now we compute all orbits of Sylp. The new element must
         permute these orbits.
      */
      allorbs(lporb, orno);
      bnoorno = orno[gbase[bno]];
      bt = 0;
      opno = bnoorno;
      opct = 1;
      skct = 0;
      for (i = 1; i <= *lporb; i++) {
        orbp[i] = 0;
        deft[i] = 0;
      }
      seek = 1;
      while (seek) {
        skct++;
        if (skct > par4) {
          printf("Exiting to try normrun.\n");
          for (i = 1; i <= npt; i++)
            if (svhptr[bno][i] == -2)
              svhptr[bno][i] = 0;
          goto exit;
        }
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
              fprintf(stderr, "Premature end of search.\n");
              return (-1);
            }
            (*expcp)--;
            bt = 1;
            adno = ntbpt[adno];
            b = 1;
            if (adno == bno) {
              im = expimage(gbase[adno]);
              l = orno[im];
              if (lporb[l] > 1)
                for (i = 1; i <= npt; i++)
                  if (orno[i] == l)
                    svhptr[adno][i] = -2;
              /* Setting sv= -2 when bno=adno avoids testing more elements
                 than necessary which map the main orbit of gbase[bno] to
                 another.
              */
            }
          }
        }
        incadno = 1;
        while (incadno) {
          cb = gbase[adno];
          im = expimage(cb);
          if (adno == bno && svhptr[adno][im] != 0) {
            ok = 0;
            break;
          }
          if (bt) {
            for (i = 1; i <= *lporb; i++) {
              j = orbp[i];
              if (deft[j] >= adno) {
                orbp[i] = 0;
                deft[j] = 0;
              }
            }
            opno = bnoorno;
            opct = 1;
            l = orbp[opno];
            while (l != 0 && l != bnoorno) {
              opct++;
              opno = l;
              l = orbp[l];
            }
          }
          bt = 0;
          deforbp();
          if (ok == 0)
            break;
          /* This means either that the element does not permute the orbits,
             or that it permutes the base orbit with a cycle of p' length.
          */
          else
            tp[cb] = im;
          if (adno == lnt) {
            incadno = 0;
            for (cb = 1; cb <= npt; cb++)
              if (invbase[cb] == 0) {
                im = expimage(cb);
                deforbp();
                if (ok == 0) {
                  bt = 1;
                  break;
                }
                else
                  tp[cb] = im;
              }
          }
          else
            adno = ntfpt[adno];
        }
        if (ok) {
          ord = findord(tp, itp);
          while (ord % prime == 0)
            ord /= prime;
          if (ord > 1) {
            for (i = 1; i <= npt; i++) {
              n = i;
              for (j = 1; j <= ord; j++)
                n = tp[n];
              itp[n] = i;
            }
            invert(itp, tp);
          }
        }
        /* if ok then tp is a p-element permuting the orbits. The final test
           is to check that it normalizes Sylp.
        */
        if (ok)
          for (i = stp; i < np2 && ok; i += 2) {
            *cp = 0;
            for (j = fnt; j <= lnt; j = ntfpt[j]) {
              k = image(tp[pptr[i][itp[gbase[j]]]]);
              if (svhptr[j][k] != 0)
                addsv(k, svhptr[j]);
              else {
                ok = 0;
                bt = 1;
                break;
              }
            }
            if (ok == 0 && ord > 1)
              for (j = bno; j < lnt; j = ntfpt[j])
                tp[gbase[j]] = expimage(gbase[j]);
          }
        /* if ok then tp is the new element of Sylp */
        if (ok) {
          seek = 0;
          fndelt();
        }
      }
    }
    for (i = 1; i <= npt; i++)
      if (svhptr[bno][i] == -2)
        svhptr[bno][i] = 0;
  }
exit:
  op = fopen(inf2, "w");
  fprintf(op, "%4d%4d%4d%4d\n", npt, *pno, nb, 4);
  printbaselo(nb, gbase, lorbdef);
  printpsv(nb, pno, svhptr);
  if (bno == 0)
    fprintf(op, "%3d", prime);
  for (i = stp; i < np2; i += 2)
    fprintf(op, "%4d", facord[i]);
  fprintf(op, "\n");
  fclose(op);
  if (x <= 0 && bno > 0) {
    op = fopen(outf1, "w");
    *pno = 0;
    for (i = 0; i < stp; i += 2)
      if (pptr[i][npt + 1] >= bno) {
        (*pno)++;
        pno[*pno] = i;
      }
    for (i = 1; i < bno; i++)
      if (lorbg[i] > 1) {
        lorbg[i] = 1;
        for (j = 1; j <= npt; j++)
          svgptr[i][j] = 0;
        svgptr[i][gbase[i]] = -1;
      }
    fprintf(op, "%4d%4d%4d%4d\n", npt, *pno, nb, 3);
    printbaselo(nb, gbase, lorbg);
    printpsv(nb, pno, svgptr);
    fclose(op);
  }
  return (bno);
}

int deforbp(void)
/* This tests whether the fact that tp sends cb (=gbase[adno]) to im
   contradicts permutation of the orbits. If so then ok is set false.
*/
{
  short i, j, k, l;
  i = orno[cb];
  j = orno[im];
  k = orbp[i];
  if (k == 0) {
    if (deft[j] == 0)
      if (lporb[i] == lporb[j]) {
        orbp[i] = j;
        deft[j] = adno;
        if (i == opno) {
          l = j;
          while (l != bnoorno && l != 0) {
            opct++;
            opno = l;
            l = orbp[l];
          }
          if (l != 0 && opct % prime != 0) {
            ok = 0;
            bt = 1;
          }
        }
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
}

int ranelt(void)
/* finds random element of group */
{
  short i, k, l;
  printf("Calling ranelt.\n");
  for (i = bno; i <= nb; i = ntfpt[i]) {
    l = rand() % lorbg[i];
    k = gorb[i][l];
    if (i == bno && svhptr[bno][k] != 0) {
      i = ntbpt[i];
      continue;
    }
    if (k != gbase[i]) {
      (*expcp)++;
      if (i <= lexp) {
        expcp[*expcp] = i;
        exprep(k, i, svgptr[i]);
      }
      else
        expcp[*expcp] = (k > gbase[i]) ? start[i] + l : start[i] + l + 1;
    }
  }
  for (i = 1; i <= npt; i++) {
    tp[i] = expimage(i);
    itp[i] = 0;
  }
}

int findord(short * p, short * ip)
/* Computes order of perm and places its inverse in ip */
{
  short ord, i, j, k, l;
  ord = 1;
  for (i = 1; i <= npt; i++)
    if (ip[i] == 0) {
      j = i;
      l = 1;
      k = p[j];
      ip[j] = j;
      while (k != i) {
        l++;
        j = k;
        k = p[j];
        ip[k] = j;
      }
      ord = lcm(ord, l);
    }
  printf("ord=%d\n", ord);
  return (ord);
}

int lcm(int x, int y)
{
  short a, b, c;
  a = x;
  b = y;
  while ((c = a % b) != 0) {
    a = b;
    b = c;
  }
  return (x * y / b);
}

int fndelt(void)
{
  short i, k, l;
  (*pno)++;
  pno[*pno] = np2 + 1;
  k = lorbdef[bno];
  l = orbitsv(gbase[bno], svhptr[bno], k);
  facord[np2 + 1] = l / k;
  lorbdef[bno] = l;
  np2 += 2;
  adno = bno;
  printf("Found element in Sylp-group. Lorbdef=%d.\n", l);
  for (i = 1; i <= npt; i++)
    if (svhptr[bno][i] == -2)
      svhptr[bno][i] = 0;
}
