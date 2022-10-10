#include "defs.h"
#include "permfns.h"

extern char  wrd, nt, isbase, inf[], outf1[], outf2[], fixed[];
extern short perm[], sv[], cp[], actgen[], orb[], base[], lorb[], order[],
    pno[], *pptr[], *svptr[], mp, mb, mnpt;
extern int psp, svsp;
short      npt;
FILE *     ip, *op;

/* The data structures in this program are similar to most permutation group
   programs. Permutations are numbered 0,1,2,3,... (where 2x+1 is usually the
   inverse of 2x) and stored in the array perm. Perm no x is pointed to by
   pptr[x]. npt=no. of points. Usually pptr[i][npt+1] gives the number of the
   group in the stabilizer chain in which the perm lies. So this no is i if
   the perm fixes the first i-1 base points. nb=no of base points. The base
   and the lengths of the orbits in the stab chain are stored in base and
   lorb. Schreier vectors are stored in sv, and pointed to by svptr[i],
   i=1,..,nb. pno is a list of *pno (=pno[0]) perm nos. cp is a similar list
   of length *cp, but it is used to represent the product of the perms
   cp[1]cp[2]..cp[*cp]. This product can be evaluated by the procedure image
   in permfns.c
*/
short gpprog(void)
{
  short i, j, k, l, m, nperms, nb, jobt, np2, seek, cord, ocord, given,
      ordknown, trivrel, lsv, u, v, w, x, y, z, mxp, mnb, bno, id;
  int   quot;
  float grporder;

  if ((ip = fopen(inf, "r")) == 0) {
    fprintf(stderr, "Cannot open file %s\n", inf);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &jobt);
  if (npt > mnpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    return (-1);
  }
  if (jobt > 0) {
    fprintf(stderr, "Wrong input format.\n");
    return (-1);
  }
  quot = psp / (npt + 1);
  if (quot > mp)
    quot = mp;
  mxp = quot;
  quot = svsp / npt;
  if (quot > mb)
    quot = mb;
  mnb = quot;
  if (mnb >= npt)
    mnb = npt - 1;
  if (nb > mnb) {
    fprintf(stderr, "nb to big. Increase SVSP (or MB).\n");
    return (-1);
  }
  /* pptr[i] is the i th permutation, svptr[i] is the i the Schreier vector.
   */
  for (i = 0; i < mxp; i++)
    pptr[i] = perm + i * (npt + 1) - 1;
  for (i = 1; i <= mnb; i++)
    svptr[i] = sv + (i - 1) * npt - 1;

  /* Next we read any base points and orbit lengths */
  if (nb != 0) {
    for (i = 1; i <= npt; i++)
      orb[i] = 0;
    for (i = 1; i <= nb; i++) {
      fscanf(ip, "%hd", base + i);
      if (orb[base[i]] != 0) {
        fprintf(stderr, "Repeated base element.\n");
        return (-1);
      }
      orb[base[i]] = 1;
    }
  }
  if (jobt < 0) {
    jobt = -jobt;
    if (jobt > nb) {
      fprintf(stderr, "jobt too big.\n");
      return (-1);
    }
    seek = jobt;
    given = jobt;
    ordknown = 1;
    for (i = 1; i <= jobt; i++)
      fscanf(ip, "%hd", order + i);
  }
  else {
    seek = 0;
    given = 0;
    ordknown = 0;
  }

  /* Now we read the generating permutations */
  np2 = 2 * nperms - 2;
  for (i = 0; i <= np2; i += 2) {
    k = i / 2 + 1;
    j = readperm(pptr[i]);
    if (j == 2) {
      fprintf(stderr, "Generator no %d is not a permutation.\n", k);
      return (-1);
    }
    if (j == 1) {
      fprintf(stderr, "Generator no %d is the identity.\n", k);
      if (nt)
        return (-1);
      nperms -= 1;
      np2 -= 2;
      i -= 2;
    }
    else {
      invert(pptr[i], pptr[i + 1]);
      actgen[i] = 1;
      x = 1;
      for (z = base[x]; x <= nb && pptr[i][z] == z; z = base[x])
        x++;
      pptr[i][npt + 1] = x;
      if (x > nb) {
        if (isbase) {
          fprintf(stderr, "Given base is not a base!\n");
          return (-1);
        }
        nb++;
        for (z = 1; pptr[i][z] == z; z++) {}
        base[nb] = z;
        printf("New base element no %d is %d.\n", nb, z);
      }
      if (x == 1)
        printf("Generator no %d moves first base point.\n", k);
      else
        printf("Generator no %d fixes first %d base point(s).\n", k, x - 1);
    }
    if (nperms == 0) {
      fprintf(stderr, "Trivial group!\n");
      return (-1);
    }
  }
  fclose(ip);

  if (wrd) {
    op = fopen(outf2, "w");
    fprintf(op, "%4d\n", nperms);
  }
  bno = nb;
  for (i = 0; i <= mnb; i++)
    lorb[i] = 0;

/* Now the main algorithm begins */
loop:
  *pno = 0;
  if (ordknown)
    ocord = (bno == nb) ? 1 : order[bno + 1];
  /* We make a list of the perm nos that fix the first bno-1 base pts */
  for (i = 0; i <= np2; i += 2) {
    if (pptr[i][npt + 1] >= bno && actgen[i] <= bno) {
      (*pno)++;
      pno[*pno] = i;
    }
  }
  lorb[bno] = orbitsv(base[bno], svptr[bno], 0);
  if (ordknown) {
    cord = (bno >= given) ? ocord * lorb[bno] : lorb[bno];
    if (bno == seek && cord == order[bno]) {
      seek--;
      printf("Order is now as given for bno = %d.\n", bno);
      bno--;
      if (bno == 0)
        goto foundorder;
      goto loop;
    }
  }
  else
    cord = 0;
  /* Now we start generating the Schreier generators that fix the firat bno
     base points, test them for membership, and add them as new generators
     if necessary.
  */
  if (*pno != 0) {
    nperms++;
    np2 += 2;
    y = np2 + 1;
    if (y >= mxp) {
      fprintf(stderr, "Out of space. Increase PSP (or MP).\n");
      return (-1);
    }
    for (i = 1; i <= npt; i++)
      fixed[i] = 0;
    for (i = 1; i < bno; i++) {
      u = base[i];
      fixed[u] = 1;
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
          id = 0;
          for (l = bno; l <= nb; l++)
            fixed[base[l]] = 0;
          for (l = bno; l <= nb; l++) {
            v = base[l];
            u = image(v);
            if (svptr[l][u] == 0)
              goto newgen;
            addsv(u, svptr[l]);
            pptr[np2][v] = v;
            pptr[y][v] = v;
            fixed[v] = 1;
          }
          l = nb + 1;
          id = 1;
        newgen:
          if (isbase == 0)
            for (k = 1; k <= npt; k++) {
              if (fixed[k] == 0) {
                u = image(k);
                pptr[np2][k] = u;
                pptr[y][u] = k;
                if (id && k != u) {
                  id = 0;
                  nb++;
                  if (nb > mnb) {
                    fprintf(stderr, "nb to big. Increase SVSP (or MB).\n");
                    return (-1);
                  }
                  base[nb] = k;
                  printf("New base point no %d is %d.\n", nb, k);
                }
              }
            }
          if (id == 0) {
            pptr[np2][npt + 1] = l;
            actgen[np2] = bno + 1;
            if (isbase)
              for (k = 1; k <= npt; k++) {
                u = image(k);
                pptr[np2][k] = u;
                pptr[y][u] = k;
              }
            printf("New generator no %d, fixing first %d base pts is:\n",
                   nperms, l - 1);
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
            if (wrd) {
              for (m = 0; m <= *cp; m++)
                fprintf(op, "%4d", cp[m]);
              fprintf(op, "\n");
            }
            bno = l;
            goto loop;
          }
        }
      }
    }
    nperms--;
    np2 -= 2;
  }
  if (ordknown && bno > given)
    order[bno] = cord;
  bno--;

foundorder:
  if (bno == 0) {
    printf("The order of the group is:\n");
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
  else
    goto loop;
  if (wrd) {
    fprintf(op, "%d\n", -1);
    fclose(op);
  }

  /* Now we output the generating perms and Schreier vectors */
  op = fopen(outf1, "w");
  fprintf(op, "%4d%4d%4d%4d\n", npt, nperms, nb, 1);
  printbaselo(nb, base, lorb);
  *pno = 0;
  for (i = 1; i <= nperms; i++) {
    (*pno)++;
    pno[*pno] = 2 * (i - 1);
  }
  printpsv(nb, pno, svptr);
  fclose(op);
  return (0);
}
