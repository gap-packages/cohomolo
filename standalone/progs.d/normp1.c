#include "defs.h"
#include "permfns.h"

extern char cent, sym, opt, hgst, nop, nonb[], inf1[], inf2[], inf3[],
    outf1[], outf2[];
extern short mp, mexp, mb, mnpt, risp, perm[], sv[], cp[], orb[], gbase[],
    hbase[], obase[], nbase[], lorbg[], lorbn[], lorbh[], ntno[], reg[],
    endorno[], ntorno[], tsv1[], tsv2[], tsv3[], genorb[], expcp[], fp[],
    pno[], start[], space[], ipno[], *pptr[], *svgptr[], *svhptr[], *svnptr[],
    *intorb[], *horno[], *hlorb[], *expptr[], *imorno[], *imlorb[],
    *orbperm[], *deftime[];
extern int psp, sp, svsp;

short npt, npt1, nptd2, nbg, nbh, stg, sth, nf, nrego, nnth, lnth, nhorbs,
    ad1, *reginfo, hgno, lexp, *norsh, *norsim, *defim, *orlist, *imhno,
    *inim, fsth, *adpt, *bindor, cb, bpt, nint, lrego, *cted, *rego, *gorno1,
    *glorb1, *gorno2, *glorb2, *test;
int  reginfolen;
char bdone, gdone;

/* normrun is a long and complex program. In this first half, we choose the
   bases to be used for H (=inf2) and G (=inf1), and record information about
   the orbit structure of H. We wish to have as many tests for
   non-normalization as possible. The search itself is in the second half
   normp2.c
*/
FILE *ip, *op;

int nprg1(void)
{
  char  heqg, inreg, con1, con2, fpt;
  short i, j, k, l, m, regb, regbno, regorep, rpno, *intno, *tsv, lo, *ptr,
      mxb, mxp, mxexp, nph, npg, maxl, ct, *spptr, *opno;
  int rsp, spn, quot;
  /* If we are storing a minimal gen set for H, we read it in */
  if (hgst) {
    if ((ip = fopen(inf3, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf3);
      return (-1);
    }
    fscanf(ip, "%hd%hd%hd%hd", &npt, &hgno, &k, &l);
    if (npt > mnpt) {
      fprintf(stderr, "npt too big. Increase NPT.\n");
      return (-1);
    }
    seeknln();
    if (k != 0)
      seeknln();
    if (l != 0)
      seeknln();
    npt1 = npt + 1;
    for (i = 0; i < hgno; i++) {
      pptr[i] = perm + npt * i - 1;
      readvec(pptr[i], 0);
      seeknln();
    }
    fclose(ip);
  }
  else
    hgno = 0;
  sth = hgno;
  /* Now read H in */
  if ((ip = fopen(inf2, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf2);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nph, &nbh, &l);
  if (npt > mnpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    return (-1);
  }
  if (l <= 2) {
    fprintf(stderr, "Wrong input format for H.\n");
    return (-1);
  }
  npt1 = npt + 1;
  spn = psp - hgno * npt;
  quot = spn / npt1;
  if (quot + hgno > mp)
    quot = mp - hgno;
  mxp = 2 * (quot / 2);
  if (2 * (nph + 1) + hgno > mxp) {
    fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
    return (-1);
  }
  ptr = perm + hgno * npt - 1;
  for (i = 0; i < mxp; i += 2) {
    pptr[i + hgno] = ptr + i * npt1;
    pptr[i + hgno + 1] = ptr + (i + 1) * npt1;
    fp[i + hgno] = i + hgno + 2;
  }
  fp[hgno + mxp - 2] = -1;
  /* The perm linker fp is required by the base changing program */
  if (sym)
    quot = svsp / npt;
  else
    quot = svsp / (2 * npt);
  if (quot > mb)
    quot = mb;
  mxb = quot;
  mb = mxb;
  if (nbh > mb) {
    fprintf(stderr, "nbh too big. Increase SVSP (or MB).\n");
    return (-1);
  }
  if (sym == 0)
    for (i = 1; i <= mxb; i++)
      svgptr[i] = sv + (i - 1) * npt - 1;
  for (i = 1; i <= mxb; i++)
    svhptr[i] = sv + svsp - i * npt - 1;
  readbaselo(nbh, hbase, lorbh);
  readpsv(sth, nbh, nph, svhptr);
  fp[sth + 2 * (nph - 1)] = -1;
  nf = sth + 2 * nph;
  fclose(ip);
  /* H is read. If sym false, read G in. */

  if (sym) {
    for (i = 1; i < npt; i++)
      lorbg[i] = npt + 1 - i;
    nbg = 0;
  }
  else {
    if ((ip = fopen(inf1, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf1);
      return (-1);
    }
    fscanf(ip, "%hd%hd%hd%hd", &k, &npg, &nbg, &l);
    stg = nf;
    if (k != npt) {
      fprintf(stderr, "npt for G and H do not agree.\n");
      return (-1);
    }
    if (nbg > mb) {
      fprintf(stderr, "nbg too big. Increase SVSP (or MB).\n");
      return (-1);
    }
    if (2 * (nph + npg + 1) + hgno > mxp) {
      fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
      return (-1);
    }
    if (l <= 2) {
      fprintf(stderr, "Wrong input format for G.\n");
      return (-1);
    }
    readbaselo(nbg, gbase, lorbg);
    readpsv(stg, nbg, npg, svgptr);
    fclose(ip);
    /* G is read. */
    fp[stg + 2 * (npg - 1)] = -1;
    nf = stg + 2 * npg;
    /* We record the original base for G as obase. gbase is variable */
    for (i = 1; i <= nbg; i++)
      obase[i] = gbase[i];
    obase[0] = nbg;
    while (lorbg[nbg] == 1)
      nbg--;
  }

  /* Now we define arrays within array sp.  reginfo records information about
     regular orbits of H, which is used to limit possibilities of
     the automorphisms that normalizing elements of H in G could induce.
     Max amount of space for this is risp.
     Choice of base elements in G and H proceeds roughly as follows.
     Set intno=1.
     a) Choose a base point for G and H in longest K-orbit O.
        (This is not always the best choice, and it can be overridden by
         using the -o option.)
        K is initially equal to H, and later is the stabilizer in H of the
        first few base points of H already chosen.
        Record the points in O in intorb[intno]. This will be used also in
        then main search, since normalizing elements must permute orbits.
     b) Choice of base points is restricted to O, until no more
        base points for G in this set are possible.
        The first stage is to choose any base points of G (but not H) that are
        possible within the fixed point set F of the  stabilizer in K of a
     point in O. (This may even be the whole of O.) This is done, in the
     program following label L3 with inreg set false. The action of the
     appropriate subgroup K of H on this F is regular, and information about
     this action is recorded in reginfo as mentioned above. When this is done,
     increment intno, and return to a), where we have now replaced K by its
     stabilizer of a point in O. In the program, we jump to label L2. c) When
     no more base points are possible inside O, we use endorno to point from
     the last such point to the number of the first. In the program, we now
     jump to L3, with inreg=true. The idea is to look for further orbits O1 of
     K, for which the kernel of the action is at least as large as that on O
     itself. If we find such orbits, then we choose as many base points of G
     (but not H of course) in O1. This saves a lot of time in the main search,
     since the automorphism of the possible normalizing element of H induced
     on the action of K on O1 can be deduced completely from that induced on
     the action of K on O. Relevant information required for such computations
     is again encoded in reginfo.

        When we are seeking the centralizer, rather than normalizer, then the
        automorphism is trivial on a zero-length orbit O, and all orbits are
        chosen as in c).
     d) When no more points can be chosen as in c), replace K by the
     appropriate stabilizer, decrement intno, and return to a).

     Let i be the number of the base point of G (not the point itself). Then:
     ntno[i]>0 means that gbase[i] is the ntno[i]-th H-base point
     (hbase[ntno[i]]) The first ad1-1 such points are not useful in search,
     since H and G have equal orbits up to this point. For i>=ad1, ntorno[i]>0
     and ntorno[i]-ntno[i]=ad1-1 (constant). These numbers are zero when
     seeking centralizer.

     If ntorno[i]>0 and within the corresponding K-orbit O there are G base
     points in F chosen as in b) above, then reg[i]<0, and -reg[i] is the
     number of this particular F. (Counted by nrego). In this case, there must
     be at least two such G-base points in F. reg[i+1] (the first of these) is
     also negative, and -reg[i+1] is the the size of F. For j>=2, and points
     within F, reg[i+j] points to an address in reginfo. This address contains
     the number i, and following addresses contain information needed for
     computing automorphisms, as mentioned above.

     For G-base points chosen as in c) above, we have reg[i]=0 for the first
     point in O1 (there maust be at least 2), and the point is recognized by
     putting ntorno[pt]= -x-1, where x is the number of the first G-base point
     in the principal orbit O. For j>=1 and points within O1, we have
     reg[i+j]>0, and the address reg[i+j] in reginfo contains -pt (negative,
     to distinguish from points as in last para), where pt is gbase[i]. Again,
     the following addresses contain information needed to compute
     automorphisms. endorno[i]>0 means i ends orbit begun at gbase no
     endorno[i].
  */
  reginfo = space;
  nhorbs = 0;
  gorno1 = risp + space;
  glorb1 = gorno1 + npt1;
  for (i = 0; i < mb; i++)
    intorb[i] = space + sp - (i + 1) * npt1;
  ptr = space + sp - (mb + 1) * npt1;
  glorb2 = ptr;
  ptr -= npt1;
  gorno2 = ptr;
  ptr -= npt;
  intno = ptr;
  ptr -= npt;
  tsv = ptr;
  ptr -= npt;
  cted = ptr;
  ptr -= npt;
  rego = ptr;
  ptr -= mp / 2;
  opno = ptr;
  rsp = ptr - space - risp - 2 * npt1;
  if (rsp <= 0) {
    fprintf(stderr, "Out of general space. Increase SPACE.\n");
    return (-1);
  }

  printf("Initial analysis and choice of base.\n");

  fpt = 0;
  nint = 0;
  gdone = 0;
  cb = 1;
  nnth = 1;
  heqg = 1;
  nrego = 0;
  nhorbs = 0;
  reginfolen = 0;
  regorep = 0;
  lnth = 0;
  intorb[0][0] = npt;
  for (i = 1; i <= npt; i++)
    cted[i] = 0;
  for (i = 1; i <= nbg; i++) {
    intorb[0][i] = gbase[i];
    cted[gbase[i]] = 1;
  }
  j = nbg;
  for (i = 1; i <= npt; i++)
    if (cted[i] == 0) {
      j++;
      intorb[0][j] = i;
    }
  for (i = 1; i < npt; i++) {
    reg[i] = 0;
    nonb[i] = 1;
  }
  nonb[npt] = 1;
  if (cent) {
    *pno = 0;
    permnos(sth, npt, nnth);
    allorbs(glorb1, gorno1);
    *opno = *pno;
    for (i = 1; i <= *pno; i++)
      opno[i] = pno[i];
    nhorbs++;
    horno[nhorbs] = gorno1;
    hlorb[nhorbs] = glorb1;
    gorno1 = glorb1 + *glorb1;
    glorb1 = gorno1 + npt1;
    rsp -= (*glorb1 + npt1);
    if (rsp <= 0) {
      fprintf(stderr, "Out of general space. Increase SPACE.\n");
      return (-1);
    }
    for (i = 1; i <= npt; i++) {
      gorno1[i] = i;
      glorb1[i] = 1;
    }
    glorb1[0] = npt;
    allorbs(glorb2, gorno2);
    skfaithorb();
    ct = 0;
    inreg = 1;
    regbno = 0;
    goto L3;
  }
  ct = 0;
  inreg = 0;
  regbno = 0;

L1:
  permnos(sth, npt, nnth);
  allorbs(glorb1, gorno1);
  if (nint > 0) {
    fpt = 1;
    ct = 0;
    lrego = 0;
    goto L3;
  }
L2:
  maxl = 1;
  if (opt) {
    ct = 0;
    for (i = 1; i <= npt; i++)
      orb[i] = 0;
    for (i = 1; i <= *intorb[nint]; i++) {
      l = intorb[nint][i];
      m = glorb1[gorno1[l]];
      if (orb[m] == 0 && m > 1) {
        ct++;
        orb[m] = 1;
        if (ct == 1) {
          bpt = l;
          maxl = m;
        }
        else {
          if (ct == 2) {
            printf("Choose one of following base points:\n");
            printf("%3d:   orbit length=%3d        ", bpt, maxl);
          }
          printf("%3d:   orbit length=%3d        ", l, m);
          if (ct % 2 == 0)
            printf("\n");
        }
      }
    }
    if (ct > 1) {
      printf("\n?      ");
      scanf("%hd", &bpt);
      maxl = glorb1[gorno1[bpt]];
    }
  }
  else
    for (i = 1; i <= *intorb[nint]; i++) {
      j = intorb[nint][i];
      k = glorb1[gorno1[j]];
      if (k > maxl) {
        bpt = j;
        maxl = k;
      }
    }
  if (maxl > 1) {
    if (newbasept(1) == -1)
      return (-1);
    nint++;
    intno[nint] = cb;
    ntno[cb] = nnth;
    k = gorno1[bpt];
    for (i = 1; i <= npt; i++)
      cted[i] = 0;
    bdone = 0;
    ct = 0;
    i = 1;
    while (i) {
      j = (bdone) ? i : gbase[i];
      if (gorno1[j] == k && cted[j] == 0) {
        cted[j] = 1;
        ct++;
        intorb[nint][ct] = j;
      }
      if (bdone)
        i = (i == npt) ? 0 : i + 1;
      else if (i == nbg) {
        bdone = 1;
        i = 1;
      }
      else
        i++;
    }
    *intorb[nint] = ct;
    lnth = cb;
    nnth++;
    if (heqg) {
      if (lorbg[cb] == lorbh[cb]) {
        if (sym == 0) {
          permnos(stg, npt, cb);
          allorbs(glorb2, gorno2);
        }
        if (sym || (*glorb2 == *glorb1)) {
          ntorno[cb] = 0;
          endorno[cb] = 0;
          if (gdone) {
            fprintf(stderr, "H = G.\n");
            return (1);
          }
          cb++;
          goto L1;
        }
      }
      else {
        ad1 = cb;
        heqg = 0;
      }
    }
    nhorbs++;
    ntorno[cb] = nhorbs;
    endorno[cb] = 0;
    horno[nhorbs] = gorno1;
    hlorb[nhorbs] = glorb1;
    gorno1 = glorb1 + *glorb1;
    glorb1 = gorno1 + npt1;
    rsp -= (*glorb1 + npt1);
    if (rsp <= 0) {
      fprintf(stderr, "Out of general space. Increase SPACE.\n");
      return (-1);
    }
    if (gdone)
      goto EXIT;
    cb++;
    goto L1;
  }

L3:
  k = (inreg) ? lrego : *intorb[nint];
  for (j = 1; j <= k; j++) {
    bpt = (inreg) ? rego[j] : intorb[nint][j];
    con1 = 1;
    i = stg;
    if (fpt) {
      if (glorb1[gorno1[bpt]] == 1)
        lrego++;
      else
        con1 = 0;
    }
    if (con1)
      con1 = (sym) ? nonb[bpt] && (gdone == 0) : pptr[i][npt1] >= cb;
    while (con1) {
      con2 = (sym) ? 1 : pptr[i][bpt] != bpt;
      if (con2) {
        if (newbasept(0) == -1)
          return (-1);
        ntno[cb] = 0;
        ntorno[cb] = 0;
        endorno[cb] = 0;
        if (heqg) {
          heqg = 0;
          ad1 = cb;
        }
        if (inreg) {
          if (ct == 0) {
            regorep = bpt;
            ntorno[cb] = -regbno - 1;
            endorno[cb] = 0;
            reg[cb] = 0;
          }
          else {
            if (ct == 1) {
              *pno = *opno;
              for (i = 1; i <= *pno; i++)
                pno[i] = opno[i];
              orbitsv(regorep, tsv, 0);
            }
            *cp = 0;
            addsv(bpt, tsv);
            reginfolen++;
            reg[cb] = reginfolen;
            reginfo[reginfolen] = -regorep;
            reginfolen++;
            reginfo[reginfolen] = *cp;
            for (i = 1; i <= *cp; i++) {
              reginfolen++;
              reginfo[reginfolen] = cp[i];
            }
            if (reginfolen > risp) {
              fprintf(stderr, "Out of reginfo space. Increase RISP.\n");
              return (-1);
            }
          }
        }
        cb++;
        ct++;
        con1 = 0;
      }
      else {
        i = fp[i];
        con1 = pptr[i][npt1] >= cb;
      }
    }
  }
  if (inreg) {
    if (gdone)
      goto EXIT;
    if (skfaithorb()) {
      ct = 0;
      goto L3;
    }
    inreg = 0;
    goto L2;
  }
  if (nint > 0) {
    regbno = intno[nint];
    if (fpt) {
      if (ct <= 1)
        reg[regbno] = 0;
      else {
        regb = gbase[regbno];
        nrego++;
        reg[regbno] = -nrego;
        reg[regbno + 1] = -lrego;
        lo = 0;
        *pno = 0;
        rpno = nf;
        for (i = 1; i <= ct; i++) {
          k = gbase[i + regbno];
          if (i == 1 || tsv[k] == 0) {
            *cp = 0;
            addsv(k, svhptr[ntno[regbno]]);
            (*pno)++;
            rpno = fp[rpno];
            pno[*pno] = rpno;
            ipno[rpno + 1] = gbase[i + regbno];
            for (j = 1; j <= npt; j++) {
              l = image(j);
              pptr[rpno][j] = l;
              pptr[rpno + 1][l] = j;
            }
            lo = orbitsv(regb, tsv, lo);
          }
          else {
            *cp = 0;
            addsv(k, tsv);
            reginfolen++;
            reg[i + regbno] = reginfolen;
            reginfo[reginfolen] = regbno;
            reginfolen++;
            reginfo[reginfolen] = *cp;
            for (j = *cp; j >= 1; j--) {
              reginfolen++;
              reginfo[reginfolen] = ipno[cp[j]];
            }
            if (reginfolen > risp) {
              fprintf(stderr, "Out of reginfo space. Increase RISP.\n");
              return (-1);
            }
          }
        }
      }
      fpt = 0;
      goto L2;
    }
    if (gdone)
      goto EXIT;
    nint--;
    if (lorbh[ntno[regbno]] >= 3) {
      if (ntorno[cb - 1] < 0)
        ntorno[cb - 1] = 0;
      endorno[cb - 1] = regbno;
      permnos(sth, nnth - 1, ntno[regbno]);
      *opno = *pno;
      for (i = 1; i <= *pno; i++)
        opno[i] = pno[i];
      allorbs(glorb2, gorno2);
      glorb2[gorno2[hbase[nnth - 1]]] = 0;
      if (skfaithorb()) {
        inreg = 1;
        ct = 0;
        goto L3;
      }
    }
    goto L2;
  }

EXIT:
  /* This concludes the initial analysis and choice of base. Now we redefine
     arrays in sp in preparation for the main search.
  */
  /* printf("nbg,nbh,lnth,ad1,nhorbs=%3d%3d%3d%3d%3d\n",
          nbg,nbh,lnth,ad1,nhorbs);
  printf("gbase:"); for (i=1;i<=nbg;i++) printf("%3d",gbase[i]);printf("\n");
  printf("hbase:"); for (i=1;i<=nbh;i++) printf("%3d",hbase[i]);printf("\n");
  printf("ntno:"); for (i=1;i<=nbg;i++) printf("%3d",ntno[i]);printf("\n");
  printf("ntorno:"); for (i=1;i<=nbg;i++)
  printf("%3d",ntorno[i]);printf("\n"); printf("endorno:"); for
  (i=1;i<=nbg;i++) printf("%3d",endorno[i]);printf("\n"); printf("reg:"); for
  (i=1;i<=nbg;i++) printf("%3d",reg[i]);printf("\n");
  */

  /* Now we start to reassign space in sp. reginfo is kept, and horno
     and hlorb are shifted back. Generators of H are copied, since they may
     be changed later, when base changes are made.
     imorno and imlorb will be the images and their lengths of the H-orbits
     under a possible normalizing element.
  */
  spptr = space + reginfolen;
  ptr = space + risp;
  while (ptr != gorno1)
    *(++spptr) = *(++ptr);
  k = risp - reginfolen;
  for (i = 1; i <= nhorbs; i++) {
    horno[i] -= k;
    hlorb[i] -= k;
  }
  test = spptr + 1;
  test[0] = npt;
  for (i = 1; i <= npt; i++)
    test[i] = i;
  spptr += npt1;
  for (i = 1; i <= nhorbs; i++) {
    k = *hlorb[i];
    imorno[i] = spptr;
    spptr += npt1;
    imlorb[i] = spptr;
    spptr += (k + 1);
    orbperm[i] = spptr;
    spptr += k;
    deftime[i] = spptr;
    spptr += k;
    for (j = 1; j <= npt; j++)
      imorno[i][j] = horno[i][j];
    for (j = 1; j <= k; j++) {
      imlorb[i][j] = hlorb[i][j];
      deftime[i][j] = 0;
      orbperm[i][j] = 0;
    }
    imlorb[i][0] = hlorb[i][0];
  }
  /* for (i=1;i<=nrego;i++) {regsv[i]=spptr; spptr +=npt;} */
  defim = spptr;
  spptr += npt;
  if (sym == 0) {
    adpt = spptr;
    spptr += nbg;
  }
  orlist = spptr;
  spptr += npt;
  inim = spptr;
  spptr += npt;
  bindor = spptr;
  spptr += npt;
  nptd2 = npt / 2;
  norsh = spptr;
  spptr += nptd2;
  norsim = spptr;
  spptr += nptd2 + 1;
  imhno = spptr;
  spptr += (mxp + 1);
  if (spptr - space >= sp) {
    fprintf(stderr, "Out of general space. Increase SPACE.\n");
    return (-1);
  }

  fsth = nf;
  for (i = sth; i != -1; i = fp[i]) {
    if (nf == -1) {
      fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
      return (-1);
    }
    imhno[i + 1] = nf + 1;
    for (j = 1; j <= npt1; j++) {
      pptr[nf][j] = pptr[i][j];
      pptr[nf + 1][j] = pptr[i][j];
    }
    k = nf;
    nf = fp[nf];
  }
  fp[k] = -1;
  for (i = 1; i <= nbg; i++)
    if (reg[i] > 0) {
      j = reg[i];
      if (reginfo[j] < 0) {
        j++;
        ptr = reginfo + j;
        k = *ptr;
        for (l = 1; l <= k; l++) {
          ptr++;
          *ptr = imhno[*ptr];
        }
      }
    }

  /* printf("reginfo:"); for (i=0;i<= *reginfo;i++) printf("%3d",reginfo[i]);
  printf("\n"); */
  /* Now we store as many cosetreps as we have room for, for use in main
   * search.
   */
  if (sym == 0) {
    spn = (sp - (spptr - space) - 1);
    printf("Remaining space for expansion = %d.\n", spn);
    mxexp = spn / npt;
    if (mxexp > mexp)
      mxexp = mexp;
    lexp = 0;
    ct = 0;
    for (i = nbg; i > ad1; i--) {
      ct += (lorbg[i] - 1);
      if (ct + i > mxexp) {
        lexp = i;
        break;
      }
    }
    if (lexp == 0)
      lexp = ad1;
    for (i = ad1; i <= lexp; i++) {
      expptr[i] = spptr;
      spptr += npt;
      start[i] = i - 1;
    }
    printf("lexp = %d.  Expanding G.\n", lexp);
    start[lexp + 1] = lexp;
    k = lexp;
    for (i = lexp + 1; i <= nbg; i++) {
      l = lorbg[i];
      start[i + 1] = start[i] + l - 1;
      if (l > 1)
        for (j = 1; j <= npt; j++)
          if (svgptr[i][j] > 0) {
            k++;
            expptr[k] = spptr;
            spptr += npt;
            exprep(j, k, svgptr[i]);
          }
    }
  }
  return (0);
}

int permnos(int no, int uf, int lf)
/* This finds the perms fixing between uf-1 and lf-1 base points, and lists
   their numbers in pno. As we are assuming format>=3, these will always be
   in order.
*/
{
  short i;
  *pno = 0;
  while (no != -1) {
    i = pptr[no][npt1];
    if (i < lf)
      no = -1;
    else if (i <= uf) {
      (*pno)++;
      pno[*pno] = no;
      no = fp[no];
    }
    else
      no = fp[no];
  }
  return (0);
}

int newbasept(int hfl)
/* Changes gbase[nnth] to bpt, If hfl=1 same for hbase[nnth] */
{
  if (hfl) {
    if (hbase[nnth] != bpt) {
      printf("H base no %d changed from %d to %d.\n", nnth, hbase[nnth], bpt);
      if ((intbase(bpt, nnth, &sth, &nbh, hbase, lorbh, svhptr)) == -1)
        return (-1);
    }
    else
      printf("H base no %d remains %d.\n", nnth, hbase[nnth]);
  }
  if (sym == 0) {
    if (gbase[cb] != bpt) {
      printf("G base no %d changed from %d to %d.\n", cb, gbase[cb], bpt);
      if ((intbase(bpt, cb, &stg, &nbg, gbase, lorbg, svgptr)) == -1)
        return (-1);
    }
    else
      printf("G base no %d remains %d.\n", cb, bpt);
  }
  else {
    printf("G base no %d defined as %d.\n", cb, bpt);
    nbg++;
    gbase[cb] = bpt;
  }
  gdone = (sym) ? cb == npt - 1 : cb == nbg;
  nonb[bpt] = 0;
  return (0);
}

int skfaithorb(void)
/* We seek a possible orbit O1, as described in a comment above. */
{
  short i, j, *k, l, n, ct;
  char  fnd;
  k = intorb[nint];
  fnd = 0;
  for (j = 1; j <= *k; j++) {
    lrego = glorb2[gorno2[*(k + j)]];
    if (lrego > 1) {
      l = gorno2[*(k + j)];
      glorb2[l] = 0;
      ct = 0;
      fnd = 1;
      for (i = 1; i <= npt; i++)
        cted[i] = 0;
      bdone = (nbg == 0);
      i = 1;
      while (i) {
        n = (bdone) ? i : gbase[i];
        if (gorno2[n] == l && cted[n] == 0) {
          ct++;
          cted[n] = 1;
          rego[ct] = n;
          if (glorb1[gorno1[n]] != 1) {
            fnd = 0;
            break;
          }
        }
        if (bdone)
          i = (i == npt) ? 0 : i + 1;
        else if (i == nbg) {
          bdone = 1;
          i = 1;
        }
        else
          i++;
      }
      if (fnd)
        return (1);
    }
  }
  return (0);
}
