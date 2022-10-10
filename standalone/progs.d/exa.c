#include "defs.h"
#include "permfns.h"

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

extern char  mult, inf1[], inf2[], inf3[], outf[];
extern short rwd[], ***scoeff[], **cdpsp[], *cpsp[], csp[], wd1[], wd2[],
    wd3[], wd4[];
extern int   cspace, psp, space, cptrsp, svsp;
extern short perm[], sv[], orb[], imsp[], *ptsp[], **simcos[], pinv[], gno[],
    base[], lorb[], pno[], *pptr[], *svptr[], cord[], invg[], rno[], cp[],
    mspace[], *vec[], **mat[], mp, cdptrsp, mwl2, mpt, mb, mdim, mpr, ptrsp,
    mm, msp, mv, mwdl, mlwdl;
int   rsp;
short endr, maxcos, *fpt, *bpt, nelim, bno, lo, **imcos, *spst, **pspst, dim,
    prime, ccos, lastd, cind, nfree, stcr, endcr, fcos, bcos, lcl, npt, np2,
    *rel, *spv, **spm, ng, rwl, nb;
short *cst, *cend, ***coeff, **cpst, ***cdpst, fullsc, clsd, lkah, conch,
    **cco, *ocst, **def;
FILE *ip, *op;

int extpprog(void)
{
  short i, j, k, l, m, n, b, l1, nperms, npt1, mxp, maxv, maxm, pos, *p, *q,
      *r;
  int     quot;
  short **dp, **dq, **dr, *pc;
  float   a;
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &k);
  if (npt > mpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    return (-1);
  }
  if (nb > mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (nb * npt > svsp) {
    fprintf(stderr, "Not enough sv space. Increase SVSP.\n");
    return (-1);
  }
  if (k <= 2) {
    fprintf(stderr, "Wrong input format.\n");
    return (-1);
  }
  npt1 = npt + 1;
  np2 = 2 * nperms - 1;
  quot = psp / npt1;
  if (quot > mp)
    quot = mp;
  mxp = quot;
  if (np2 >= mxp) {
    fprintf(stderr, "Too many permutations. Increase MP and/or PSP.\n");
    return (-1);
  }
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + i * npt1 - 1;
  for (i = 1; i <= nb; i++)
    svptr[i] = sv + (i - 1) * npt - 1;
  readbaselo(nb, base, lorb);
  readpsv(0, nb, nperms, svptr);
  fclose(ip);
  for (i = 0; i < np2; i += 2) {
    invg[i] = i + 1;
    invg[i + 1] = i;
  }

  /* We need to read in matrices for generators when action is nontrivial */
  if (mult == 0) {
    if ((ip = fopen(inf3, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf3);
      return (-1);
    }
    fscanf(ip, "%hd%hd%hd", &prime, &dim, &k);
    if (prime > mpr) {
      fprintf(stderr, "Prime too big. Increase MPR.\n");
      return (-1);
    }
    if (dim > mdim) {
      fprintf(stderr, "Dim too big. Increase MDIM.\n");
      return (-1);
    }
    if (k != nperms) {
      fprintf(stderr, "No of mats wrong.\n");
      return (-1);
    }
    setpinv();
    maxv = msp / dim;
    if (maxv > mv)
      maxv = mv;
    for (i = 0; i < maxv; i++)
      vec[i] = mspace - 1 + i * dim;
    maxm = maxv / dim - 1;
    if (maxm >= mm)
      maxm = mm - 1;
    if (np2 >= maxm) {
      fprintf(stderr, "Not enough mat space. Increase MSP and/or MV,MM.\n");
      return (-1);
    }
    for (i = 0; i <= maxm; i++)
      mat[i] = vec - 1 + i * dim;
    spm = mat[maxm];
    spv = spm[1];
    for (i = 0; i < np2; i += 2) {
      readmat(mat[i]);
      inv(mat[i], mat[i + 1]);
    }
    fclose(ip);
  }

  printf("Input a,b,l1.\nmaxcos = a*lorb+b\n");
  printf("Exit lookahead after eliminating abs(l1) cosets.\n");
  printf("l1 negative for consistency checking.\n");
  scanf("%f%hd%hd", &a, &b, &l1);
  if (l1 < 0) {
    l1 = -l1;
    conch = 1;
  }
  else
    conch = 0;

  /* Now we read in the relations */
  if ((ip = fopen(inf2, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf2);
    return (-1);
  }
  fscanf(ip, "%hd%hd", &i, &j);
  if (i != nb || (mult == 0 && j != dim)) {
    fprintf(stderr, "Error in %s, line 1.\n", inf2);
    return (-1);
  }
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%hd", rno + i);
  if (mult) {
    dim = j;
    for (i = 1; i <= dim; i++)
      fscanf(ip, "%hd", cord + i);
  }
  rel = imsp;
  spst = rel;
  for (i = 1; i <= 2 * rno[1]; i++) {
    fscanf(ip, "%hd", spst);
    l = *spst;
    for (j = 1; j <= l; j++)
      fscanf(ip, "%hd", spst + j);
    spst += (l + 1);
  }
  fclose(ip);

  /* Now we are ready to start the modified Todd-Coxeter procedure. This
     follows the ordinary tcrun with procedures added in the appropriate
     places to compute the coefficients as words in the generators of G and M.
     The Tood-Coxeter procedures are in exb.c and the word operation
     procedures in exc.c
  */
  pos = 0;
  endr = -1;
  pspst = ptsp;
  cst = csp;
  ocst = cst;
  cend = csp + cspace - 1;
  cpst = cpsp;
  cdpst = cdpsp;
  rsp = cspace;
  for (bno = nb; bno >= 1; bno--) {
    printf("bno=%d.\n", bno);
    l = rno[bno] - pos;
    for (i = 1; i <= 2 * l; i++) {
      endr++;
      endr += rel[endr];
    }
    pos = rno[bno];
    lo = lorb[bno];
    if (lo == 1)
      continue;
    maxcos = lo * a + b;
    ng = 0;
    tmalloc(def, short *, maxcos + 1);
    for (i = 0; i < np2; i += 2)
      if (pptr[i][npt1] >= bno)
        ng += 2;
      else
        break;
    k = (ng + 1) * maxcos;
    l = k + maxcos + 1;
    k -= lo;
    if (cpst + k - cpsp > cptrsp) {
      fprintf(stderr, "Out of cptrsp. Increase CPTRSP.\n");
      return (-1);
    }
    if (spst + l - imsp > space) {
      fprintf(stderr, "Out of im space. Increase SPACE.\n");
      return (-1);
    }
    if (pspst + ng - ptsp > ptrsp) {
      fprintf(stderr, "Out of ptrsp. Increase PTRSP.\n");
      return (-1);
    }
    if (cdpst + ng - cdpsp > cdptrsp) {
      fprintf(stderr, "Out of cdptr Increase CDPTRSP.\n");
      return (-1);
    }
    imcos = pspst;
    pspst += ng;
    gno[bno] = ng;
    for (i = 0; i < ng; i++) {
      imcos[i] = spst - 1;
      p = imcos[i];
      spst += maxcos;
      while (++p < spst)
        *p = 0;
    }
    fpt = spst;
    bpt = spst + maxcos;
    coeff = cdpst;
    cdpst += ng;
    for (i = 0; i < ng; i++) {
      coeff[i] = cpst - 1;
      dp = cpst - 1;
      cpst += maxcos;
      while (++dp < cpst)
        *dp = 0;
    }
    cco = cpst - lo - 1;
    *pno = 0;
    for (i = 0; i < np2; i += 2)
      if (pptr[i][npt1] >= bno) {
        (*pno)++;
        pno[*pno] = i;
        if (pptr[i][npt1] > bno) {
          imcos[i][1] = 1;
          imcos[i + 1][1] = 1;
          coeff[i][1] = cst;
          *cst = 1;
          *(++cst) = i;
          *(++cst) = 0;
          cst++;
          coeff[i + 1][1] = cst;
          *cst = 1;
          *(++cst) = i + 1;
          *(++cst) = 0;
          cst++;
          rsp -= 6;
        }
      }
    if (orbitsv(base[bno], svptr[bno], 0) != lo) {
      fprintf(stderr, "Orbit length wrong!\n");
      return (-1);
    }
    if (mult == 0) {
      tmalloc(def[1], short, 1);
      def[1][0] = 0;
    }
    for (l = 1; l <= lo; l++) {
      n = orb[l];
      cp[n] = l;
      if ((k = svptr[bno][n]) != -1) {
        j = pptr[k][n];
        m = cp[j];
        imcos[k - 1][m] = l;
        imcos[k][l] = m;
        if (mult == 0) {
          tmalloc(def[l], short, def[m][0] + 2);
          def[l][0] = def[m][0] + 1;
          def[l][1] = k;
          for (i = 1; i <= def[m][0]; i++)
            def[l][i + 1] = def[m][i];
        }
      }
    }

    lastd = lo;
    cind = lo;
    nfree = lo + 1;
    for (i = 1; i <= lo; i++)
      bpt[i] = i - 1;
    for (i = 0; i < maxcos; i++)
      fpt[i] = i + 1;
    fpt[lo] = 0;
    fpt[maxcos] = 0;
    ccos = 1;
    lkah = 0;
    /* At last we are ready to proceed with the enumeration */
    while (ccos != 0) {
      clsd = 1;
      endcr = -1;
      while (endcr != endr) {
        if (scanrel() == -1)
          return (-1);
        if (fullsc == 0) {
          clsd = 0;
          if (lkah == 0) {
            lcl = bpt[ccos];
            lkah = 1;
            nelim = 0;
            printf("Entering lookahead.\n");
          }
        }
      }
      ccos = fpt[ccos];
      if (lkah) {
        if (nelim >= l1 || ccos == 0) {
          if (cind != maxcos) {
            printf("Exiting lookahead. cind=%d.\n", cind);
            ccos = fpt[lcl];
            lkah = 0;
          }
          else {
            fprintf(stderr, "maxcos too small.\n");
            return (-1);
          }
        }
      }
    } /* while ccos!=0 */
      /* End of enumeration. Tidy up and prepare for next bno */
    p = imcos[0] + lo;
    q = imcos[1];
    k = maxcos - lo;
    dp = coeff[0] + lo;
    dq = coeff[1];
    for (i = 1; i < ng; i++) {
      imcos[i] = p;
      r = p + lo;
      while (++p <= r)
        *p = *(++q);
      q += k;
      coeff[i] = dp;
      dr = dp + lo;
      while (++dp <= dr)
        *dp = *(++dq);
      dq += k;
    }
    spst = p + 1;
    cpst = dp + 1;
    simcos[bno] = imcos;
    scoeff[bno] = coeff;
    if (bno > 1) {
      maxcos = lo;
      gc();
      ocst = cst;
    }
    if (mult == 0)
      for (i = 1; i <= lo; i++)
        tfree(def[i]);
    tfree(def);
  } /* for (bno=;... */

  /* End of program. Now output */
  op = fopen(outf, "w");
  fprintf(op, "%4d%4d\n", nb, dim);
  for (i = 1; i <= nb; i++)
    fprintf(op, " %3d", lorb[i]);
  fprintf(op, "\n");
  if (mult) {
    for (i = 1; i <= dim; i++)
      fprintf(op, "%4d", cord[i]);
    fprintf(op, "\n");
  }
  for (i = 1; i <= nb; i++)
    fprintf(op, "%4d", gno[i]);
  fprintf(op, "\n");
  for (i = 1; i <= nb; i++)
    if ((l = lorb[i]) > 1)
      for (j = 0; j < gno[i]; j++) {
        p = simcos[i][j];
        for (k = 1; k <= l; k++)
          fprintf(op, " %3d", p[k]);
        fprintf(op, "\n");
        dp = scoeff[i][j];
        for (k = 1; k <= l; k++) {
          if ((pc = dp[k]) == 0)
            fprintf(op, "%3d %3d ", 0, 0);
          else {
            fprintf(op, "%3d ", *pc);
            for (m = 1; m <= *pc; m++)
              fprintf(op, " %2d", pc[m]);
            pc += (1 + *pc);
            fprintf(op, " %3d ", *pc);
            for (m = 1; m <= *pc; m++)
              fprintf(op, " %2d", pc[m]);
          }
        }
        fprintf(op, "\n");
      }
  fclose(op);
  return (0);
}
