#include "defs.h"

extern char slg, con, check, inf1[], inf2[], inf3[], inf4[], inf5[], inf6[],
    outf1[];
extern int   msp, psp, svsp;
extern short mm, mv, mp, mpt, mpr, mdim, mb, mspace[], *vec[], **mat[],
    pinv[], perm[], sv[], cp[], *pptr[], *svptr[], base[], lorb[];
short prime, dim, *spv, **spm, npt, nb, maxv, maxm, maxp, nmat, nm2;
FILE *ip, *op;

void seeknln(void)
{
  while (getc(ip) != '\n')
    ;
}

/* We repeat a few procedures from permfns.c, since that file is not loaded
   in matcalc.
*/
void addsv(int pt, short * sv)
{
  short pn;
  pn = sv[pt];
  while (pn != -1) {
    (*cp)++;
    cp[*cp] = pn;
    pt = pptr[pn][pt];
    pn = sv[pt];
  }
}

void invert(short * ptr1, short * ptr2)
{
  short i;
  for (i = 1; i <= npt; i++)
    ptr2[ptr1[i]] = i;
}

int image(int pt)
{
  short i;
  for (i = 1; i <= *cp; i++)
    pt = pptr[cp[i]][pt];
  return (pt);
}

void readvec(short * ptr, int e)
{
  short i;
  for (i = 1; i <= npt + e; i++)
    fscanf(ip, "%hd", ptr + i);
}

void readbaselo(int nb, short * base, short * lorb)
{
  short i;
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%hd", base + i);
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%hd", lorb + i);
}

void readpsv(int e, int nb, int nperms, short ** svptr)
{
  short i, j, *k;
  for (i = 1; i <= nperms; i++) {
    j = e + 2 * i - 2;
    readvec(pptr[j], 1);
    invert(pptr[j], pptr[j + 1]);
  }
  for (i = 1; i <= nb; i++) {
    readvec(svptr[i], 0);
    for (j = 1; j <= npt; j++) {
      k = svptr[i] + j;
      if (*k > 0)
        *k += e;
    }
  }
}

void setpinv(void)
{
  short i, j;
  for (i = 0; i < prime; i++)
    pinv[i] = 0;
  for (i = 1; i < prime; i++)
    if (pinv[i] == 0)
      for (j = 1; j < prime; j++)
        if (i * j % prime == 1) {
          pinv[i] = j;
          pinv[j] = i;
          break;
        }
}

int mcprog(void)
{
  short i, j, k, l, n, nperms, np2, *p, **q;
  int   quot;
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd", &prime, &dim, &nmat);
  if (prime > mpr) {
    fprintf(stderr, "prime too big. Increase MPR.\n");
    return (-1);
  }
  if (dim > mdim) {
    fprintf(stderr, "dim too big. Increase MDIM.\n");
    return (-1);
  }
  setpinv();
  /* Set pointers to define vectors and matrices. Matrix elements are stored
     in mspace, as vectors of vectors.
  */
  quot = msp / dim;
  if (quot > mv)
    quot = mv;
  maxv = quot;
  for (i = 0; i < maxv; i++)
    vec[i] = mspace - 1 + i * dim;
  maxm = maxv / dim - 1;
  if (maxm >= mm)
    maxm = mm - 1;
  for (i = 0; i <= maxm; i++)
    mat[i] = vec - 1 + i * dim;
  /* spm is a spare matrix used for inverting, etc. similarly spv */
  spm = mat[maxm];
  spv = spm[1];
  nm2 = nmat * 2 - 1;
  if (nm2 >= maxm) {
    fprintf(stderr, "Too many mats. Increase MM.\n");
    return (-1);
  }
  for (i = 0; i < nm2; i += 2) {
    readmat(mat[i]);
    if (inv(mat[i], mat[i + 1]) == -1)
      return (-1);
  }
  fclose(ip);
  if ((ip = fopen(inf2, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf2);
    return (-1);
  }
  fscanf(ip, "%hd", &i);
  if (i != nmat) {
    fprintf(stderr, "nmat does not agree in %s.\n", inf2);
    return (-1);
  }
  /* Compute matrices for strong generators of G. Then save only those in
     gpname.sg (unless -o set)
  */
  while (1) {
    fscanf(ip, "%hd", cp);
    if (*cp == -1)
      break;
    for (i = 1; i <= *cp; i++)
      fscanf(ip, "%hd", cp + i);
    nmat++;
    nm2 += 2;
    if (nm2 >= maxm) {
      fprintf(stderr, "Too many mats. Increase MM.\n");
      return (-1);
    }
    prod(cp, mat[nm2 - 1]);
    inv(mat[nm2 - 1], mat[nm2]);
  }
  fclose(ip);
  if ((ip = fopen(inf3, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf3);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &l);
  if (npt > mpt) {
    fprintf(stderr, "npt too big. Increase MPT.\n");
    return (-1);
  }
  if (nb >= mb) {
    fprintf(stderr, "nb too big. Increase MB.\n");
    return (-1);
  }
  if (nb * npt > svsp) {
    fprintf(stderr, "svsp too small. Increase SVSP.\n");
    return (-1);
  }
  np2 = 2 * nperms - 1;
  if (l <= 0 || (l != 3 && slg) || np2 > nm2 || (slg == 0 && np2 != nm2)) {
    fprintf(stderr, "%s has invalid data nm2,np2,l=%d,%d,%d.\n", inf3, nm2,
            np2, l);
    return (-1);
  }
  quot = psp / npt;
  if (quot > mp)
    quot = mp;
  maxp = quot;
  if (np2 >= maxp) {
    fprintf(stderr, "Need more perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  for (i = 0; i < maxp; i++)
    pptr[i] = perm + i * (npt + 1) - 1;
  for (i = 1; i <= nb; i++)
    svptr[i] = sv + (i - 1) * npt - 1;
  readbaselo(nb, base, lorb);
  readpsv(0, nb, nperms, svptr);
  if (slg) {
    p = mspace + msp - 1 - nmat;
    for (i = 1; i <= nmat; i++)
      p[i] = 0;
    for (i = 1; i <= nperms; i++) {
      fscanf(ip, "%hd", cp + i);
      p[cp[i]] = 1;
    }
    fclose(ip);
    k = 0;
    for (i = 1; i <= nperms; i++) {
      q = vec - 1 + (2 * cp[i] - 2) * dim;
      mat[k] = q;
      mat[k + 1] = q + dim;
      k += 2;
    }
    for (i = 1; i <= nmat; i++)
      if (p[i] == 0) {
        q = vec - 1 + (2 * i - 2) * dim;
        mat[k] = q;
        mat[k + 1] = q + dim;
        k += 2;
      }
    nmat = nperms;
    nm2 = np2;
    /* If -t set, we check that the matrices satisfy the relations of G */
    if (check) {
      if ((ip = fopen(inf6, "r")) == 0) {
        fprintf(stderr, "Cannot open %s.\n", inf6);
        return (-1);
      }
      fscanf(ip, "%hd%hd", &i, &n);
      seeknln();
      for (i = 1; i <= n; i++) {
        fscanf(ip, "%hd", cp);
        for (j = 1; j <= *cp; j++)
          fscanf(ip, "%hd", cp + j);
        q = mat[nm2 + 1];
        prod(cp, q);
        for (k = 1; k <= dim; k++)
          for (l = 1; l <= dim; l++)
            if ((l == k && q[k][l] != 1) || (l != k && q[k][l] != 0)) {
              fprintf(stderr, "Matrices do not satisfy group relations.\n");
              fprintf(stderr, "   (Relation %d in %s.)\n", i, inf6);
              fclose(ip);
              return (-1);
            }
      }
      fclose(ip);
    }
  }
  return (0);
}

int conprog(int con)
{
  short i, j, k, n, y, fac, lm3, bn, lmat, nperms, class, *d1, *d2, *wt, *p,
      *q, *r, *cv, *cve, *cc, *sw, *bcf, *ibcf, **cbm, **icbm, **newmat;
  int  sum;
  char id;
  if (con == 1) {
    for (i = 0; i < nm2; i += 2) {
      trans(mat[i], mat[i + 1]);
      inv(mat[i + 1], mat[i]);
    }
  }
  if ((ip = fopen(inf4, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf4);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &i, &nperms, &j, &k);
  seeknln();
  if (i != npt) {
    fprintf(stderr, "npt error in %s.\n", inf4);
    return (-1);
  }
  y = con == 1 ? nperms : 1;
  if (nm2 + y >= maxp) {
    fprintf(stderr, "Need more perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  y = con == 1 ? nperms + 4 : 1;
  if (nm2 + y >= maxm) {
    fprintf(stderr, "Need more mat space. Increase MSP (or MM or MV).\n");
    return (-1);
  }
  if (j > 0)
    seeknln();
  if (k != 0)
    seeknln();
  /* In case con=1, we reverse the order of the generators, since
     the generators of P as permutation group are in reverse order from
     PCP generators.
     We deal with the case con=1 separately, since we need to read in all
     of the generators from inf4, and compute their matrices, in order to
     compute the necessary basis changes.
  */
  if (con == 1) {
    for (i = nperms; i >= 1; i--) {
      readvec(pptr[nm2 + i], 0);
      seeknln();
    }
    fclose(ip);

    /* Now we express the perms in terms of the generators of G, and compute
       their matrices.
    */
    for (i = 1; i <= nperms; i++) {
      p = pptr[nm2 + i];
      *cp = 0;
      for (j = 1; j <= nb; j++)
        addsv(image(p[base[j]]), svptr[j]);
      q = cp + *cp;
      p = q + 1;
      *p = *cp;
      while (q > cp) {
        *(++p) = (*q % 2 == 0) ? *q + 1 : *q - 1;
        q--;
      }
      p = cp + *cp + 1;
      prod(p, mat[nm2 + i]);
    }
    lmat = nm2 + nperms;
    /* Now do the basis changes.
       The base change matrices will be output to inf6, for use later in
       nqrun.
    */
    op = fopen(inf6, "w");
    fprintf(op, "%4d%4d%4d\n", prime, dim, 2);
    cbm = mat[lmat + 1];
    if (dim == 1) {
      cbm[1][1] = id = 1;
      icbm = NULL;
      bcf = NULL;
      ibcf = NULL;
      lm3 = 0;
    }
    else {
      icbm = mat[lmat + 2];
      bcf = mat[lmat + 4][1];
      lm3 = lmat + 3;
      for (i = 1; i <= dim; i++)
        bcf[i] = 0;
      ibcf = mat[lmat + 4][2];
      id = 1;
      for (n = dim; n >= 1; n--) {
        cv = mat[lm3][1];
        cve = cv + dim;
        cc = mat[lm3][2];
        for (i = dim; i >= 1; i--)
          if (bcf[i] == 0) {
            if (n == 1) {
              p = cbm[1];
              q = p + dim;
              while (++p <= q)
                *p = 0;
              cbm[1][i] = 1;
              break;
            }
            else {
              bn = i;
              p = cv;
              while (++p <= cve)
                *p = 0;
              cv[i] = 1;
              while (1) {
                for (j = 1; j <= nperms; j++) {
                  comm(cv, cc, mat[nm2 + j]);
                  for (k = dim; k > n; k--)
                    if ((fac = cc[ibcf[k]]) != 0) {
                      p = cc;
                      r = p + dim;
                      q = cbm[k] + 1;
                      while (++p <= r) {
                        sum = *p - (fac * *q);
                        *p = sum % prime;
                        if (*p < 0)
                          *p += prime;
                        q++;
                      }
                    }
                  for (k = dim; k > 0; k--)
                    if (cc[k] != 0) {
                      bn = k;
                      id = 0;
                      sw = cv;
                      cv = cc;
                      cc = sw;
                      cve = cv + dim;
                      break;
                    }
                  if (k > 0)
                    break;
                }
                if (j > nperms)
                  break;
              }
              bcf[bn] = n;
              ibcf[n] = bn;
              p = cbm[n];
              q = p + dim;
              r = cv + 1;
              fac = pinv[cv[bn]];
              while (++p <= q) {
                sum = *r * fac;
                *p = sum % prime;
                r++;
              }
              break;
            }
          }
      }
    }
    printmat(cbm);
    if (id)
      printf("No first basis change.\n");
    else { /* printf("First dual basis change matrix (new in terms of
           old):\n"); for (i=1;i<=dim;i++) {for (j=1;j<=dim;j++)
           printf("%3d",cbm[i][j]); printf("\n"); } */
      inv(cbm, icbm);
      newmat = mat[lm3];
      *cp = 3;
      cp[1] = lmat + 1;
      cp[3] = lmat + 2;
      for (i = 0; i <= lmat; i++) {
        cp[2] = i;
        prod(cp, newmat);
        mat[lm3] = mat[i];
        mat[i] = newmat;
        newmat = mat[lm3];
      }
    }
    if (dim == 1) {
      d1 = mat[lmat + 2][1];
      d2 = mat[lmat + 3][1];
      wt = mat[lmat + 4][1];
      wt[1] = class = 1;
      d1[1] = d2[1] = 0;
    }
    else {
      d1 = bcf;
      d2 = ibcf;
      wt = (dim == 2) ? mat[lmat + 5][1] : mat[lmat + 4][3];
      id = cbdef(nm2 + 1, nm2 + nperms, lmat + 1, d1, d2, wt, &class);
    }
    printmat(cbm);
    fclose(op);
    if (id)
      printf("No second basis change.\n");
    else { /* printf("Second dual basis change matrix (new in terms of
           old):\n"); for (i=1;i<=dim;i++) {for (j=1;j<=dim;j++)
           printf("%3d",cbm[i][j]); printf("\n"); } */
      inv(cbm, icbm);
      newmat = mat[lm3];
      *cp = 3;
      cp[1] = lmat + 1;
      cp[3] = lmat + 2;
      for (i = 0; i <= lmat; i++) {
        cp[2] = i;
        prod(cp, newmat);
        mat[lm3] = mat[i];
        mat[i] = newmat;
        newmat = mat[lm3];
      }
    }
    op = fopen(outf1, "w");
    fprintf(op, "%4d%4d%4d%3d\n", prime, dim, nperms, class);
    for (i = 1; i <= dim; i++)
      fprintf(op, "%3d ", wt[i]);
    fprintf(op, "\n");
    for (i = 1; i <= dim; i++)
      fprintf(op, "%3d ", d1[i]);
    fprintf(op, "\n");
    for (i = 1; i <= dim; i++)
      fprintf(op, "%3d ", d2[i]);
    fprintf(op, "\n");
    for (i = 1; i <= nperms; i++)
      printmat(mat[nm2 + i]);
    fclose(op);
  } /* con=1 */
  else {
    op = fopen(outf1, "w");
    fprintf(op, "%4d%4d%4d\n", prime, dim, nperms);
    for (i = 1; i <= nperms; i++) {
      p = pptr[nm2 + 1];
      readvec(p, 0);
      *cp = 0;
      seeknln();
      for (j = 1; j <= nb; j++)
        addsv(image(p[base[j]]), svptr[j]);
      if (con == 0) {
        q = cp + *cp;
        p = q + 1;
        *p = *cp;
        while (q > cp) {
          *(++p) = (*q % 2 == 0) ? *q + 1 : *q - 1;
          q--;
        }
        p = cp + *cp + 1;
      }
      else
        p = cp;
      prod(p, mat[nm2 + 1]);
      printmat(mat[nm2 + 1]);
    }
    fclose(ip);
    fclose(op);
  }
  return (0);
}
