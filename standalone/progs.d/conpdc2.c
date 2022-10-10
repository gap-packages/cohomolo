#include "defs.h"

extern char inf0[], inf1[], inf2[], inf3[], outf1[], outf2[], outf3[],
    temp1[], temp2[], ofstr[], ngn[], expg, exph, cr, dcr, hg, triv, cong;
extern int trsp, trptr;
extern int nfuse, tree[], perm[], gorb[], lorbg[], lorbh[], base[], fpt[],
    bpt[], coh_index[], *cp[], *trad[], *sv[], *tailad[], **cpad[],
    **svgptr[], **svhptr[], npt, **lcp, **ucp, *stop, *pst, *pend, nb, fnt,
    lnt, ind, coind, oind, cind, bno, pt, **expp, npg;
extern FILE *ip, *op;

/* The basic aim of this part of the program is to compute permutations of G
   as permutations on the left cosets of H. This is used to compute double
   coset reps if required.
*/
int cnprg2(void)
{
  int i, j, k, l, m, n, ncp, pno, mexh, lorb, ndcr, nntb, *orb, rep, *q, **p,
      **conp, *len, *orbptr, *cperm, flct, inch;
  char ngth, fuse;
  int  c;
  if (exph) {
    mexh = (pend - pst) / npt;
    lcp = cp;
    j = 0;
    for (i = 1; i <= lnt; i++)
      j += (lorbh[i] - 1);
    if (j < mexh)
      for (i = 1; i <= lnt; i++) {
        for (j = 1; j <= npt; j++)
          expp[j] = 0;
        expp[base[i]] = stop;
        p = svhptr[i];
        for (j = 1; j <= npt; j++)
          if (p[j] != 0 && p[j] != stop) {
            expp[j] = pst;
            ucp = lcp - 1;
            addsvf(j, p);
            for (k = 1; k <= npt; k++)
              pst[k] = image(k);
            pst += npt;
          }
        for (j = 1; j <= npt; j++)
          p[j] = expp[j];
      }
    else {
      printf("Not enough room to expand H.\n");
      exph = 0;
    }
  }
  conp = cp + 2 * npt;
  if (cong) {
    strcat(inf1, ".cp");
    op = fopen(inf1, "w");
    fprintf(op, "%4d%4d\n", ind, npg);
    printf("Computing coeffs of perms in G.\n");
    for (pno = 1; pno <= npg; pno++) {
      *conp = perm + (pno - 1) * 2 * npt - 1;
      printf("Coeffs of permutation no %d.\n", pno);
      /* The basic conversion algorithm now follows. In this case, coefficents
         are being output, as well as the perms themselves.
      */
      ucp = conp;
      q = tree;
      for (bno = fnt; bno != -1; bno = fpt[bno]) {
        tailad[bno] = q;
        q += 2;
        cpad[bno] = ucp + 1;
      }
      for (cind = 1; cind <= ind; cind++) {
        if (cind > 1)
          advance();
        lcp = conp;
        trptr = 0;
        for (bno = 1; bno <= nb; bno++) {
          if (fpt[bno] != 0) {
            i = tree[trptr];
            j = image(i);
            while (svhptr[bno][j] == 0) {
              trptr = tree[trptr + 1];
              i = tree[trptr];
              j = image(i);
            }
            if (expg) {
              q = svgptr[bno][i];
              if (q != stop) {
                lcp--;
                *lcp = q + npt;
              }
            }
            else
              addsvb(i, svgptr[bno]);
            trptr += 2;
          }
          else
            j = image(base[bno]);
          if (exph) {
            q = svhptr[bno][j];
            if (q != stop) {
              ucp++;
              *ucp = q;
            }
          }
          else
            addsvf(j, svhptr[bno]);
          fprintf(op, " %d", j);
        }
        fprintf(op, " %6d  ", tree[trptr]);
      } /* for (cind=1;... */
      fprintf(op, "\n");
    } /* for (pno=1;...  */
    fclose(op);
  } /* if (cong) */
  if (hg == 0)
    return (0);

  *conp = expg ? perm - 1 : pst - 1;
  flct = -1;
  printf("Converting permutations from %s.\n", inf3);
  if ((ip = fopen(inf3, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf3);
    return (-1);
  }
  fscanf(ip, "%d%d%d%d", &i, &ncp, &j, &k);
  if (k != 0)
    k = 3;
  else if (j > 0)
    k = 2;
  else
    k = 1;
  if (i != npt) {
    fprintf(stderr, "npt wrong in %s.\n", inf3);
    return (-1);
  }
  for (i = 1; i <= k; i++)
    while ((c = getc(ip)) != '\n')
      ;
  op = fopen(temp1, "w");
  fprintf(op, "%4d%4d%4d%4d\n", ind, ncp, 0, 0);

  ngth = 0;
  while (1) {
    for (pno = 1; pno <= ncp; pno++) {
      printf("Converting perm no %d.\n", pno);
      q = *conp + 1;
      for (i = 1; i <= npt; i++) {
        fscanf(ip, "%d", q);
        q++;
      }
      while ((c = getc(ip)) != '\n')
        ;
      ucp = conp;
      q = tree;
      for (bno = fnt; bno != -1; bno = fpt[bno]) {
        tailad[bno] = q;
        q += 2;
        cpad[bno] = ucp + 1;
      }
      for (cind = 1; cind <= ind; cind++) {
        if (cind > 1)
          advance();
        lcp = conp;
        trptr = 0;
        for (bno = 1; bno <= lnt; bno++) {
          if (fpt[bno] != 0) {
            i = tree[trptr];
            j = image(i);
            while (svhptr[bno][j] == 0) {
              trptr = tree[trptr + 1];
              i = tree[trptr];
              j = image(i);
            }
            if (expg) {
              q = svgptr[bno][i];
              if (q != stop) {
                lcp--;
                *lcp = q + npt;
              }
            }
            else
              addsvb(i, svgptr[bno]);
            trptr += 2;
          }
          else
            j = image(base[bno]);
          if (exph) {
            q = svhptr[bno][j];
            if (q != stop) {
              ucp++;
              *ucp = q;
            }
          }
          else
            addsvf(j, svhptr[bno]);
        }
        fprintf(op, " %3d", tree[trptr]);
      }
      fprintf(op, "\n");
    } /* for (pno=1;... */
    fclose(op);
    fclose(ip);
    flct++;
    if (flct > nfuse)
      break;
    if (flct == 0) {
      strcpy(inf2, inf0);
      strcat(inf2, ngn);
      if ((ip = fopen(inf2, "r")) == 0) {
        flct++;
        if (flct > nfuse)
          break;
        ngth = 0;
      }
      else
        ngth = 1;
    }
    if (flct > 0) {
      strcpy(inf2, inf0);
      strcat(inf2, "dcr");
      ofstr[0] = flct + '0';
      strcat(inf2, ofstr);
    }
    strcpy(outf3, inf2);
    strcat(outf3, ".nr");
    if ((ip = fopen(inf2, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", inf2);
      return (-1);
    }
    fscanf(ip, "%d%d", &i, &ncp);
    while (getc(ip) != '\n')
      ;
    if (i != npt) {
      fprintf(stderr, "npt wrong in %s.\n", inf2);
      return (-1);
    }
    op = fopen(outf3, "w");
    fprintf(op, "%3d\n", ncp);
    printf("Converting permutations from %s.\n", inf2);
  } /* while(1) */

  /* We are now ready to compute orbits of the new perms if necessary, in
     order to compute double coset reps. The tree information is no longer
     required, and the tree space is used to store the new permutations.
  */
  if (dcr) {
    if (trsp < 4 * ind) {
      fprintf(stderr, "No room to compute orbits. Increase TRSP.\n");
      return (-1);
    }
    ndcr = ind - 1;
    orb = tree - 1;
    len = orb + ind;
    orbptr = len + ind;
    cperm = orbptr + ind;
    for (i = 1; i <= ind; i++) {
      orb[i] = i;
      len[i] = 1;
      orbptr[i] = 0;
    }
    printf("Computing orbits of con perms from %s.\n", outf1);
    ip = fopen(temp1, "r");
    fscanf(ip, "%d%d", &i, &ncp);
    while ((c = getc(ip)) != '\n') {}
    flct = -1;
    fuse = 0;
    while (1) {
      for (i = 1; i <= ncp; i++) {
        q = cperm;
        for (j = 1; j <= ind; j++) {
          q++;
          fscanf(ip, "%d", q);
        }
        for (j = 1; j <= ind; j++) {
          k = orb[j];
          l = orb[cperm[j]];
          if (k != l) {
            if (fuse) {
              if (len[l] < len[k]) {
                m = k;
                k = l;
                l = m;
              }
            }
            else
              len[k] += len[l];
            len[l] = 0;
            n = k;
            while ((m = orbptr[n]) != 0)
              n = m;
            orbptr[n] = l;
            while (l != 0) {
              orb[l] = k;
              l = orbptr[l];
            }
            ndcr--;
          }
        }
      }
      fclose(ip);
      flct++;
      if (flct == 0)
        unlink(temp1);
      else
        unlink(outf3);
      if (flct > nfuse)
        break;
      if (flct == 0) {
        fuse = 1;
        if (ngth == 0)
          flct++;
        if (flct > nfuse)
          break;
      }
      strcpy(outf3, inf0);
      if (flct == 0)
        strcat(outf3, ngn);
      else {
        strcat(outf3, "dcr");
        ofstr[0] = flct + '0';
        strcat(outf3, ofstr);
      }
      strcat(outf3, ".nr");
      ip = fopen(outf3, "r");
      fscanf(ip, "%d", &ncp);
      printf("Computing orbits of con perms from %s.\n", outf3);
    }
    printf("Computing dcreps.\n");
    if (len[1] != 1)
      printf("Warning. Orb of 1 is greater than one.\n");
    nntb = 0;
    for (i = fnt; i != -1; i = fpt[i])
      nntb++;
    lcp = cp;
    ip = fopen(temp2, "r");
    op = fopen(outf2, "w");
    fprintf(op, "%4d %3d%4d%4d\n", npt, ndcr, 0, 0);
    for (rep = 2; rep <= ind; rep++)
      if ((lorb = len[rep]) > 0) {
        ucp = lcp - 1;
        for (bno = fnt; bno != -1; bno = fpt[bno]) {
          fscanf(ip, "%d", &inch);
          if (expg) {
            q = svgptr[bno][inch];
            if (q != stop) {
              ucp++;
              *ucp = q;
            }
          }
          else
            addsvf(inch, svgptr[bno]);
        }
        if (npt >= 1000)
          for (i = 1; i <= npt; i++)
            fprintf(op, "%5d", image(i));
        else
          for (i = 1; i <= npt; i++)
            fprintf(op, "%4d", image(i));
        fprintf(op, "   %4d\n", lorb);
      }
      else
        for (i = 1; i <= nntb; i++)
          fscanf(ip, "%d", &inch);
    fclose(op);
    fclose(ip);
    unlink(temp2);
  }
  /* If we are keeping the converted permutations, then invert them, to
    get action on right cosets
  */
  if (hg && dcr == 0) {
    ip = fopen(temp1, "r");
    op = fopen(outf1, "w");
    fscanf(ip, "%d%d", &i, &ncp);
    while ((c = getc(ip)) != '\n')
      ;
    if (triv)
      fprintf(op, "%4d%4d%4d%4d\n1\n", ind, ncp, 1, 0);
    else
      fprintf(op, "%4d%4d%4d%4d\n", ind, ncp, 0, 0);
    for (i = 1; i <= ncp; i++) {
      for (j = 1; j <= ind; j++) {
        fscanf(ip, "%d", &k);
        tree[k] = j;
      }
      for (j = 1; j <= ind; j++)
        fprintf(op, " %3d", tree[j]);
      fprintf(op, "\n");
    }
    fclose(ip);
    fclose(op);
    unlink(temp1);
  }
  return (0);
}
