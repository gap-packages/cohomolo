#include "defs.h"
#include "permfns.h"

#define SCP 100
extern char heqg, words, fullg, check, outf[], inf1[], inf2[], inf3[],
    outf0[];
extern short perm[], sv[], cp[], actgen[], orb[], base[], lorbc[], lorbg[],
    lorbh[], order[], pno[], *pptr[], *svptr1[], *svptr2[], *svptr3[], scpf[],
    scps[], sadpt[], mp, mb;
extern int psp, svsp;
short      npt, mxp, sth, stcom, np2, nb, npt1, cps, cpf;
char       hing, cthere, hsvth;
FILE *     ip, *op;

int optprog(void)
{
  short i, j, k, nperms, jobt, *ptr, *iptr;
  int   quot;
  char  c, err, opt[80];

  /* This program stores up to three groups G,H and C. G must always have a
     s.g. set and Schreier Vectors svptr1. H and C may have Schreier vectors
     svptr2 and svptr3. C is not always defined.
     heqg=1 means H=G. cthere=1 means C is defined.
     Perms in G have nos. from 0 to sth-2.
     Perms in H have nos. from sth to stcomm-2 (unless H=G, when H is not
     stored) stcomm-2=np2 if C not defined. If C defined, perms in C have nos
     from stcomm to np2. First read inf1 = G
  */
  if ((ip = fopen(inf1, "r")) == 0) {
    printf("Cannot open %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &jobt);
  if (jobt <= 0) {
    fprintf(stderr, "Wrong input format for G.\n");
    return (-1);
  }
  npt1 = npt + 1;
  quot = psp / npt1;
  if (quot > mp)
    quot = mp;
  mxp = quot;
  np2 = 2 * nperms - 2;
  if (spacecheck() == -1)
    return (-1);
  if (svsp / npt < 3 * nb) {
    fprintf(stderr, "Not enough sv space.\n");
    return (-1);
  }
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + i * (npt + 1) - 1;
  for (i = 1; i <= nb; i++)
    svptr1[i] = sv + (i - 1) * npt - 1;
  for (i = 1; i <= nb; i++)
    svptr2[i] = sv + (i + nb - 1) * npt - 1;
  for (i = 1; i <= nb; i++)
    svptr3[i] = sv + (i + 2 * nb - 1) * npt - 1;
  readbaselo(nb, base, lorbg);
  readpsv(0, nb, nperms, svptr1);
  sth = (heqg) ? 0 : 2 * nperms;
  fclose(ip);
  /* if fullg read from gpname.inperm to determine minimal generating set for
   * G */
  if (fullg) {
    if ((ip = fopen(inf3, "r")) == 0) {
      printf("Cannot open %s.\n", inf3);
      return (-1);
    }
    fscanf(ip, "%hd%hd", &i, &k);
    fclose(ip);
  }
  else
    k = nperms;
  /* actgen=1 means we have a genuine generator of G */
  for (i = 1; i <= k; i++)
    actgen[2 * i - 2] = 1;
  for (i = k + 1; i <= nperms; i++)
    actgen[2 * i - 2] = 2;
  if (heqg == 0)
  /* Read inf2 = H */
  {
    if ((ip = fopen(inf2, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf2);
      return (-1);
    }
    fscanf(ip, "%hd%hd%hd%hd", &i, &nperms, &j, &k);
    err = 0;
    np2 += 2 * nperms;
    if (spacecheck() == -1)
      return (-1);
    if (i != npt || (k > 0 && j != nb))
      err = 1;
    if (err == 0)
      if (k > 0) {
        for (i = 1; i <= nb; i++) {
          fscanf(ip, "%hd", &j);
          if (j != base[i]) {
            err = 1;
            break;
          }
        }
        if (err == 0) {
          for (i = 1; i <= nb; i++)
            fscanf(ip, "%hd", lorbh + i);
          readpsv(sth, nb, nperms, svptr2);
          hsvth = 1;
        }
        for (i = 1; i <= nperms; i++)
          actgen[2 * i + sth - 2] = 1;
      }
      else {
        seeknln();
        if (j != 0)
          seeknln();
        if (k != 0)
          seeknln();
        for (i = 1; i <= nperms; i++) {
          j = sth + 2 * i - 2;
          ptr = pptr[j];
          iptr = ptr + npt1;
          readvec(ptr, 0);
          seeknln();
          invert(ptr, iptr);
          actgen[j] = 1;
          for (k = 1; k <= nb; k++)
            if (ptr[base[k]] != base[k])
              break;
          ptr[npt1] = k;
        }
        hsvth = 0;
      }
    if (err) {
      fprintf(stderr, "Wrong input format for H.\n");
      return (-1);
    }
    fclose(ip);

    /* Check whether H is a subgroup of G and express generators of H as words
       in gens of G if required.
    */
    hing = 1;
    for (i = 1; i <= nperms; i++) {
      *cp = 1;
      cp[1] = sth + 2 * i - 2;
      if (test(svptr1, words) == 0) {
        hing = 0;
        if (words == 0)
          break;
      }
    }
  }
  else {
    sth = 0;
    hing = 1;
    hsvth = 1;
    for (i = 1; i <= nb; i++)
      lorbh[i] = lorbg[i];
  }
  if (hing) {
    if (check)
      printf("true\n");
    else
      printf("H is a subgroup of G.\n");
  }
  else {
    if (check)
      printf("false\n");
    else
      printf("H is not a subgroup of G.\n");
  }
  if (check)
    return (0);
  cthere = 0;
  stcom = sth + 2 * nperms;
  np2 = stcom - 2;

  while (1)
  /* Now start processing according to chosen options. */
  {
    printf("Choose option. ('l' for list options.)\n");
    scanf("%s", opt);
    if (strcmp(opt, "ex") == 0)
      scanf("%hd", &k);
    while ((c = getchar()) != '\n')
      ;
    if (strcmp(opt, "l") == 0) {
      printf("Options:\nnc    Put C=normal closure of H in G,\n");
      printf("comm  Put C=[G,H],\nint   Put C=G ^ H,\njoin  Put C=<G,H>,\n");
      printf("core  Replace H by its core w.r.t. G,\n");
      printf("oph   Output H,\nopc   Output C,\n");
      printf("hsv   Calculate s.g.set, Schreier vectors and order of H,\n");
      printf("rgh   Replace G by H and set H=G,\n");
      printf("rghc  Replace G by H and H by C,\n");
      printf("rhc   Replace H by C,\n");
      printf("ex n  Test if perm. no n of H lies in G. If so express it in "
             "gens of G,\n");
      printf("q     Quit.\n");
      printf("WARNING: core,hsv,rghc and rhc render C undefined.\n");
    }
    else if (strcmp(opt, "nc") == 0) {
      if (heqg)
        printf("H = G. Normal closure = G.\n");
      else if (nc(0, 0, sth - 2) == -1)
        return (-1);
    }
    else if (strcmp(opt, "comm") == 0) {
      if (comm() == -1)
        return (-1);
    }
    else if (strcmp(opt, "int") == 0) {
      if (hing)
        printf("H <= G. Intersection = H.\n");
      else if (hsvth == 0)
        printf("hsv must be called first.\n");
      else if (intsect(svptr1, svptr2, svptr3, stcom) == -1)
        return (-1);
    }
    else if (strcmp(opt, "core") == 0) {
      if (heqg)
        printf("H = G. Core = G.\n");
      else if (hsvth == 0)
        printf("hsv must be called first.\n");
      else if (core() == -1)
        return (-1);
    }
    else if (strcmp(opt, "join") == 0) {
      if (hing)
        printf("H <= G. Join = G.\n");
      else
        join();
    }
    else if (strcmp(opt, "oph") == 0) {
      if (hsvth == 0)
        fprintf(stderr, "H has no sv. Output not done.\n");
      else if (outp(0) == -1)
        return (-1);
    }
    else if (strcmp(opt, "opc") == 0) {
      if (cthere == 0)
        fprintf(stderr, "No subgroup C to output.\n");
      else if (outp(1) == -1)
        return (-1);
    }
    else if (strcmp(opt, "hsv") == 0) {
      if (hsvth)
        fprintf(stderr, "H has sv already.\n");
      else {
        cthere = 0;
        np2 = stcom - 2;
        if (gp(sth, nb, lorbh, svptr2) == -1)
          return (-1);
        hsvth = 1;
        stcom = np2 + 2;
      }
    }
    else if (strcmp(opt, "rhc") == 0) {
      if (cthere == 0)
        fprintf(stderr, "No subgroup C.\n");
      else
        rhc();
    }
    else if (strcmp(opt, "rgh") == 0) {
      if (heqg)
        fprintf(stderr, "H = G already.\n");
      else
        rgh(0);
    }
    else if (strcmp(opt, "rghc") == 0) {
      if (heqg)
        fprintf(stderr, "H = G already.\n");
      else if (cthere == 0)
        printf("No subgroup C.\n");
      else
        rgh(1);
    }
    else if (strcmp(opt, "ex") == 0) {
      if (heqg)
        fprintf(stderr, "H=G, so this option is pointless.\n");
      else if (abs(k) > nperms || k == 0)
        fprintf(stderr, "Inappropriate perm no.\n");
      else {
        j = (k > 0) ? sth + 2 * (k - 1) : sth + 2 * (-k - 1) + 1;
        *cp = 1;
        cp[1] = j;
        if (test(svptr1, 0) == 0)
          printf("Perm is not in G.\n");
        else {
          printf("Perm is in G. It is:\n");
          for (i = *cp; i > 1; i--) {
            j = cp[i];
            k = (j % 2 == 0) ? -j / 2 - 1 : (j - 1) / 2 + 1;
            printf("%d", k);
            if (i == 2)
              printf("\n");
            else
              printf(" * ");
          }
        }
      }
    }
    else if (strcmp(opt, "q") == 0)
      break;
    else
      printf("Invalid option.\n");
  } /* Main option choosing loop  */
  return (0);
}

int gp(int spno, int sbno, short * lorb, short ** svptr)
/* Compute s.g. set and group order of a group K (=H or C).
   Similar to program gprun.
   spno=no. of first perm in K; the last no. is the external np2.
   lorb and svptr for placing orbit lengths and Schreier vectors of K.
   sbno=no. of base points.
   Warning: It is assumed that the given base really is a base for K,
            and no new base points are added. This is the main difference
            from gprun. The answer will be false if this assumption is false.
            However, the call from intsect uses a deliberately smaller value
            of sbno, to save time.

*/
{
  short i, j, l, m, bno, u, v, w, x, y, z, lsv;
  char  trivrel, id;
  float grporder;
  bno = sbno;
loop:
  *pno = 0;
  z = (spno == 0) ? sth - 2 : np2;
  for (i = spno; i <= z; i += 2) {
    if (pptr[i][npt1] >= bno && actgen[i] <= bno) {
      (*pno)++;
      pno[*pno] = i;
    }
  }
  lorb[bno] = orbitsv(base[bno], svptr[bno], 0);
  if (*pno != 0) {
    np2 += 2;
    y = np2 + 1;
    if (spacecheck() == -1)
      return (-1);
    for (i = 1; i < bno; i++) {
      u = base[i];
      pptr[np2][u] = u;
      pptr[y][u] = u;
    }
    for (i = 1; i <= lorb[bno]; i++) {
      *cp = 0;
      addsv(orb[i], svptr[bno]);
      for (w = 1, x = *cp; w <= x; w++, x--) {
        if (w == x)
          cp[w] -= 1;
        else {
          u = cp[w];
          cp[w] = cp[x] - 1;
          cp[x] = u - 1;
        }
      }
      lsv = *cp;
      for (j = 1; j <= *pno; j++) {
        *cp = lsv;
        u = pno[j] + 1;
        trivrel = (*cp > 0) ? cp[*cp] == (u - 1) : 0;
        if (trivrel == 0) {
          (*cp)++;
          cp[*cp] = u;
          id = 1;
          for (l = bno; l <= nb; l++) {
            v = base[l];
            u = image(v);
            if (svptr[l][u] == 0) {
              id = 0;
              break;
            }
            addsv(u, svptr[l]);
            pptr[np2][v] = v;
            pptr[y][v] = v;
          }
          if (id == 0) {
            pptr[np2][npt1] = l;
            actgen[np2] = bno + 1;
            for (m = 1; m <= npt; m++) {
              u = image(m);
              pptr[np2][m] = u;
              pptr[y][u] = m;
            }
            printf("New generator fixing first %d base pts is:\n", l - 1);
            for (m = 1; m <= *cp; m++) {
              z = cp[m];
              y = z / 2;
              x = (z == 2 * y) ? y + 1 : -(y + 1);
              printf("%3d", x);
              if (m == *cp)
                printf("\n");
              else
                printf("*");
            }
            bno = l;
            goto loop;
          }
        }
      }
    }
    np2 -= 2;
  }
  bno--;

  if (bno == 0) {
    if (spno != 0) {
      printf("Present order of the group is:\n");
      for (i = 1; i <= nb; i++) {
        printf("%3d", lorb[i]);
        if (i == nb)
          printf(" =\n");
        else
          printf("*");
      }
      grporder = 1;
      for (i = 1; i <= nb; i++)
        grporder *= lorb[i];
      printf("%g\n", grporder);
    }
    return (0);
  }
  else
    goto loop;
}

int test(short ** svp, int exp)
/* Test whether the current perm (in cp) lies in the group G, using the
   Schreier Vector svp. If exp is true, print out word for element in
   generators of G.
*/
{
  short i, k, l;
  for (l = 1; l <= nb; l++) {
    k = image(base[l]);
    if (svp[l][k] == 0) {
      if (exp)
        printf("Permutation is not in group.\n");
      return (0);
    }
    addsv(k, svp[l]);
  }
  if (exp)
    for (i = cp[0]; i > 1; i--) {
      k = cp[i];
      l = (k % 2 == 0) ? -(k + 2) / 2 : (k + 1) / 2;
      printf("%2d", l);
      if (i > 2)
        printf(" *");
      else
        printf("\n");
    }
  return (1);
}

int nc(int th, int sno, int eno)
/* Compute normal closure of a group K under perms from sno to eno.
   If th=0, first copy H to C, put th=stcomm and np2 the last one.
   In any case, K is then generated by perms from th to np2.
*/
{
  short i, j, l, m, n, *p, *q, *r;
  if (th == 0)
  /* Copy H and its Schreier vector into C if necessary. */
  {
    p = pptr[sth];
    q = pptr[stcom];
    r = q;
    while (++p <= r)
      *(++q) = *p;
    np2 = 2 * stcom - sth - 2;
    if (spacecheck() == -1)
      return (-1);
    p = actgen + sth - 1;
    q = actgen + stcom - 1;
    r = q;
    while (++p <= r)
      *(++q) = *p;
    if (hsvth) {
      p = svptr2[1];
      q = svptr3[1];
      r = q;
      l = stcom - sth;
      while (++p <= r)
        *(++q) = (*p <= 0) ? *p : *p + l;
      p = lorbh;
      q = lorbc;
      r = p + nb;
      while (++p <= r)
        *(++q) = *p;
    }
    else if (gp(stcom, nb, lorbc, svptr3) == -1)
      return (-1);
    th = stcom;
  }
  for (i = th; i <= np2; i += 2)
    if (actgen[i] == 1)
      for (j = sno; j <= eno; j += 2)
        if (actgen[j] == 1) {
          *cp = 3;
          cp[1] = j + 1;
          cp[2] = i;
          cp[3] = j;
          if (test(svptr3, 0) == 0)
          /* Add new gen to C */
          {
            np2 += 2;
            p = pptr[np2];
            q = pptr[np2 + 1];
            *cp = 3;
            if (spacecheck() == -1)
              return (-1);
            for (n = 1; n <= npt; n++) {
              m = image(n);
              p[n] = m;
              q[m] = n;
            }
            actgen[np2] = 1;
            for (n = 1; n <= nb; n++)
              if (p[base[n]] != base[n])
                break;
            p[npt1] = n;
            if (gp(stcom, n, lorbc, svptr3) == -1)
              return (-1);
          }
        }
  cthere = 1;
  return (0);
}

int comm(void)
/* Compute commutator subgroup C = [G,H] */
{
  short i, j, k, *p, *q, *r, m, n, gnext, hnext;
  p = svptr3[1];
  r = p + nb * npt;
  while (++p <= r)
    *p = 0;
  np2 = stcom - 2;
  if (spacecheck() == -1)
    return (-1);
  for (i = 1; i <= nb; i++) {
    svptr3[i][base[i]] = -1;
    lorbc[i] = 1;
  }
  k = heqg ? stcom : sth;
  /* First define C to be generated by commutators of generators */
  for (i = 0; i < k; i += 2)
    if (actgen[i] == 1)
      for (j = sth; j < stcom; j += 2)
        if (actgen[j] == 1) {
          *cp = 4;
          cp[1] = i + 1;
          cp[2] = j + 1;
          cp[3] = i;
          cp[4] = j;
          if (test(svptr3, 0) == 0) {
            np2 += 2;
            p = pptr[np2];
            q = pptr[np2 + 1];
            *cp = 4;
            if (spacecheck() == -1)
              return (-1);
            for (n = 1; n <= npt; n++) {
              m = image(n);
              p[n] = m;
              q[m] = n;
            }
            actgen[np2] = 1;
            for (n = 1; n <= nb; n++)
              if (p[base[n]] != base[n])
                break;
            p[npt1] = n;
            if (gp(stcom, n, lorbc, svptr3) == -1)
              return (-1);
          }
        }
  if (np2 < stcom) {
    printf("[G,H] is trivial.\n");
    cthere = 0;
    return (0);
  }
  gnext = stcom;
  hnext = stcom;
  /* Now take normal closure under G and H */
  while (1) {
    if (gnext > np2)
      break;
    n = heqg ? stcom - 2 : sth - 2;
    if (nc(gnext, 0, n) == -1)
      return (-1);
    gnext = np2 + 2;
    if (hing || hnext > np2)
      break;
    if (nc(hnext, sth, stcom - 2) == -1)
      return (-1);
  }
  return (0);
}

int outp(int c)
/* Out put a group. H if c=0, C if c=1. */
{
  short nperms, i, s, e, **sv;
  char  sn[12];
  printf("Input name of subgroup to be output:    ");
  scanf("%s", sn);
  strcpy(outf, outf0);
  strcat(outf, sn);
  op = fopen(outf, "w");
  nperms = c ? (np2 + 2 - stcom) / 2 : (stcom - sth) / 2;
  fprintf(op, "%4d%4d%4d%4d\n", npt, nperms, nb, 1);
  if (c)
    printbaselo(nb, base, lorbc);
  else
    printbaselo(nb, base, lorbh);
  s = c ? stcom : sth;
  e = c ? np2 : stcom - 2;
  *pno = 0;
  for (i = s; i <= e; i += 2)
    pno[++(*pno)] = i;
  sv = c ? svptr3 : heqg ? svptr1 : svptr2;
  printpsv(nb, pno, sv);
  fclose(op);
  return (0);
}

int intsect(short ** sv1, short ** sv2, short ** sv3, int stint)
/* Compute intersection of groups with Schreier vectors sv1 and sv2.
   This uses a becktrack search thro' the sv1 group.
   The result has S. vector sv3, and s.g. set from stint to np2.
   It is necessary for the current perm to expand in both directions, so
   we start it at position SCP defined above. Chaos would result if this
   should be too small. In general it goes from cps to cpf. in the array cp.
*/
{
  short i, m, n, k, adno, adpt, bno, *p, *q, *r;
  np2 = stint - 2;
  p = sv3[1];
  r = p + nb * npt;
  while (++p <= r)
    *p = 0;
  if (spacecheck() == -1)
    return (-1);
  for (i = 1; i <= nb; i++) {
    sv3[i][base[i]] = -1;
    lorbc[i] = 1;
  }
  scps[1] = SCP;
  scpf[1] = scps[1] - 1;
  for (bno = nb; bno >= 1; bno--) {
    adno = bno;
    cps = scps[1];
    cpf = scpf[1];
    scps[bno] = cps;
    scpf[bno] = cpf;
    adpt = 0;
    while (1)
    /* First advance thro' search */
    {
      adpt++;
      while (adpt <= npt && sv1[adno][adpt] == 0)
        adpt++;
      if (adpt > npt) {
        if (adno == bno)
          break;
        adno--;
        cps = scps[adno];
        cpf = scpf[adno];
        adpt = sadpt[adno];
        continue;
      }
      if (adno == bno && (sv2[bno][adpt] == 0 || sv3[bno][adpt] != 0))
        continue;
      addsvb(adpt, sv1[adno]);
      k = im(base[adno]);
      if (sv2[adno][k] == 0) {
        cps = scps[adno];
        continue;
      }
      addsvf(k, sv2[adno]);
      sadpt[adno] = adpt;
      if (adno == nb)
      /* We have a new generator in the intersection. */
      {
        np2 += 2;
        if (spacecheck() == -1)
          return (-1);
        p = pptr[np2];
        q = pptr[np2 + 1];
        p[npt1] = bno;
        actgen[np2] = 1;
        cps = scps[bno];
        for (n = 1; n <= npt; n++) {
          m = im(n);
          p[n] = m;
          q[m] = n;
        }
        if (gp(stint, bno, lorbc, sv3) == -1)
          return (-1);
        adno = bno;
        adpt = sadpt[bno];
        cpf = scpf[bno];
        continue;
      }
      adno++;
      scps[adno] = cps;
      scpf[adno] = cpf;
      adpt = 0;
    }
  }
  if (np2 < stint) {
    printf("Intersection is trivial.\n");
    cthere = 0;
    return (0);
  }
  else
    cthere = 1;
  return (1);
}

int im(int pt)
/* Image under current perm used in intsect */
{
  short i;
  for (i = cps; i <= cpf; i++)
    pt = pptr[cp[i]][pt];
  return (pt);
}

int addsvb(int pt, short * sv)
/* Similar to addsv in permfns.c, but goes backwards, using current perm
   in intsect
*/
{
  short pn;
  pn = sv[pt];
  while (pn != -1) {
    cps--;
    cp[cps] = pn - 1;
    pt = pptr[pn][pt];
    pn = sv[pt];
  }
}

int addsvf(int pt, short * sv)
/* Same but goes forwards */
{
  short pn;
  pn = sv[pt];
  while (pn != -1) {
    cpf++;
    cp[cpf] = pn;
    pt = pptr[pn][pt];
    pn = sv[pt];
  }
}

int core(void)
/* Replace H by its core under G */
{
  short i, j, l, m, n, stint, *p, *q, *r;
  np2 = stcom - 2;
  if (spacecheck() == -1)
    return (-1);
restart:
  for (i = 0; i < sth; i += 2)
    if (actgen[i] == 1)
      for (j = sth; j < stcom; j += 2) {
        *cp = 3;
        cp[1] = i + 1;
        cp[2] = j;
        cp[3] = i;
        if (test(svptr2, 0) == 0)
        /* The conjugate of the generator of H under the generator of G does
           not lie in H, so we must take an intersection. First compute g-1 H
           g and then its intersection with H
        */
        {
          for (j = sth; j < stcom; j += 2) {
            *cp = 3;
            cp[2] = j;
            np2 += 2;
            p = pptr[np2];
            q = pptr[np2 + 1];
            for (n = 1; n <= npt; n++) {
              m = image(n);
              p[n] = m;
              q[m] = n;
            }
            actgen[np2] = 1;
            for (n = 1; n <= nb; n++)
              if (p[base[n]] != base[n])
                break;
            p[npt1] = n;
          }
          if (gp(stcom, nb, lorbc, svptr3) == -1)
            return (-1);
          stint = np2 + 2;
          if ((n = intsect(svptr2, svptr3, svptr1, stint)) == -1)
            return (-1);
          if (n == 0)
          /* If core is trivial, we exit after putting heqg=1 for good
             measure. */
          {
            printf("Core is trivial. Setting H=G.\n");
            if (gp(0, nb, lorbg, svptr1) == -1)
              return (-1);
            stcom = sth;
            np2 = sth - 2;
            sth = 0;
            hing = 1;
            heqg = 1;
            cthere = 0;
            return (0);
          }
          l = stint - sth;
          p = svptr1[1];
          q = svptr2[1];
          r = p + nb * npt;
          while (++p <= r)
            *(++q) = (*p <= 0) ? *p : *p - l;
          for (i = 1; i <= nb; i++)
            lorbh[i] = lorbc[i];
          p = pptr[stint];
          q = pptr[sth];
          r = pptr[np2] + 2 * npt1;
          while (++p <= r)
            *(++q) = *p;
          np2 = sth + np2 - stint;
          stcom = np2 + 2;
          for (i = sth; i <= np2; i++)
            actgen[i] = 1;
          goto restart;
        }
      }
  /* We have messed up svptr1 during this computation, so we recompute it */
  cthere = 0;
  gp(0, nb, lorbg, svptr1);
  *cp = 1;
  hing = 1;
  for (i = sth; i <= np2; i += 2) {
    cp[1] = i;
    if (test(svptr1, 0) == 0) {
      hing = 0;
      break;
    }
  }
  if (hing)
    printf("H is a subgroup of G.\n");
  else
    printf("H is not a subgroup of G.\n");
  return (0);
}

int rhc(void)
/* Replace H by C */
{
  short i, l, *p, *q, *r;
  cthere = 0;
  if (heqg) {
    sth = stcom;
    heqg = 0;
  }
  l = stcom - sth;
  p = svptr3[1];
  q = svptr2[1];
  r = p + nb * npt;
  while (++p <= r)
    *(++q) = (*p <= 0) ? *p : *p - l;
  for (i = 1; i <= nb; i++)
    lorbh[i] = lorbc[i];
  if (heqg == 0) {
    p = pptr[stcom];
    q = pptr[sth];
    r = pptr[np2] + 2 * npt1;
    while (++p <= r)
      *(++q) = *p;
    p = actgen + stcom - 1;
    q = actgen + sth - 1;
    r = actgen + np2 + 1;
    while (++p <= r)
      *(++q) = *p;
  }
  np2 = sth + np2 - stcom;
  stcom = np2 + 2;
  hsvth = 1;
  *cp = 1;
  hing = 1;
  for (i = sth; i <= np2; i += 2) {
    cp[1] = i;
    if (test(svptr1, 0) == 0) {
      hing = 0;
      break;
    }
  }
  if (hing)
    printf("H is a subgroup of G.\n");
  else
    printf("H is not a subgroup of G.\n");
}

int rgh(int c)
/* Replace G by H, and is c=1 also H by C */
{
  short i, l, *p, *q, *r;
  l = sth;
  p = svptr2[1];
  q = svptr1[1];
  r = p + nb * npt;
  while (++p <= r)
    *(++q) = (*p <= 0) ? *p : *p - l;
  for (i = 1; i <= nb; i++)
    lorbg[i] = lorbh[i];
  p = pptr[sth];
  q = pptr[0];
  r = pptr[stcom];
  while (++p <= r)
    *(++q) = *p;
  p = actgen + sth - 1;
  q = actgen - 1;
  r = actgen + stcom - 1;
  while (++p <= r)
    *(++q) = *p;
  if (c) {
    sth = stcom - sth;
    rhc();
  }
  else {
    stcom -= sth;
    np2 = stcom - 2;
    sth = 0;
    cthere = 0;
    hsvth = 1;
    hing = 1;
    heqg = 1;
  }
}

int join(void)
/* Compute C=<G,H> */
{
  short i, *p, *q, *r;
  np2 = stcom - 2;
  for (i = 0; i < stcom; i += 2)
    if (actgen[i] == 1) {
      if (spacecheck() == -1)
        return (-1);
      np2 += 2;
      p = pptr[i];
      q = pptr[np2];
      r = p + 2 * npt1;
      if (spacecheck() == -1)
        return (-1);
      while (++p <= r)
        *(++q) = *p;
      actgen[np2] = 1;
    }
  gp(stcom, nb, lorbc, svptr3);
  cthere = 1;
  return (0);
}

int spacecheck(void)
{
  if (np2 + 1 >= mxp) {
    fprintf(stderr, "Out of perm space. Increase PSP (or MP).\n");
    return (-1);
  }
  return (0);
}
