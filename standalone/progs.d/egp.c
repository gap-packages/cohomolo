#include "defs.h"
#include "permfns.h"

extern char  inf[], inf2[], outf[], full, stabcall;
extern short perm[], sv2[], cp[], orb[], base[], lorb[], fixpt[], fixb[],
    pno[], *pptr[], *svptr[], mnpt, mp, mb;

short npt;
FILE *ip, *op;

int egprog(void)
{
  short nperms, ngens, nb, stab, olfb, nlfb, i, l, m, z;
  if (full) {
    if ((ip = fopen(inf2, "r")) == 0) {
      printf("Cannot open %s.\n", inf2);
      return (-1);
    }
    fscanf(ip, "%hd%hd", &i, &ngens);
    fclose(ip);
  }
  if ((ip = fopen(inf, "r")) == 0) {
    printf("Cannot open %s.\n", inf);
    return (-1);
  }
  if (stabcall) {
    printf("Stabilize how many points?    ");
    scanf("%hd", &stab);
  }
  else
    stab = 0;
  fscanf(ip, "%hd%hd%hd%hd", &npt, &nperms, &nb, &l);
  if (npt > mnpt) {
    fprintf(stderr, "npt,too big. Increase NPT.\n");
    return (-1);
  }
  if (nperms > mp) {
    fprintf(stderr, "nperms too big. Increase MXP.\n");
    return (-1);
  }
  if (nb > mb) {
    fprintf(stderr, "nb too big. Increase MNB.\n");
    return (-1);
  }
  if (l <= 0) {
    fprintf(stderr, "Wrong input format.\n");
    return (-1);
  }
  for (i = 0; i <= nperms; i++)
    pptr[i] = perm + (npt + 1) * i - 1;
  for (i = 1; i <= nb; i++)
    svptr[i] = sv2 + npt * (i - 1) - 1;
  for (i = 1; i <= npt; i++)
    pptr[0][i] = i;
  readbaselo(nb, base, lorb);
  for (i = 1; i <= nperms; i++)
    readvec(pptr[i], 1);
  fclose(ip);

  /* Main algorithm begins */
  op = fopen(outf, "w");
  *pno = 0;
  fixpt[nb + 1] = 0;
  for (l = nb; l > stab; l--) {
    for (i = 1; i <= nperms; i++)
      if (pptr[i][npt + 1] == l) {
        (*pno)++;
        pno[*pno] = i;
      }
    fixpt[l] = *pno;
  }
  for (l = stab; l >= 1; l--)
    fixpt[l] = fixpt[stab + 1];
  for (l = nb; l > stab; l--) {
    olfb = fixpt[l + 1];
    nlfb = fixpt[l];
    if ((nlfb - olfb) > 1)
      for (i = olfb + 1; i <= nlfb; i++) {
        z = pno[i];
        /* We test whether perm no z is redundant as a generator.
           If so, we leave pno[i]=0, if not we put it back to z.
        */
        if (full == 0 || z > ngens) {
          pno[i] = 0;
          for (m = l; m > stab && pno[i] == 0; m--) {
            *pno = fixpt[m];
            if (orbitsv(base[m], svptr[m], 0) < lorb[m])
              pno[i] = z;
          }
        }
      }
  }
  *pno = 0;
  l = 0;
  for (i = nb; i >= 1; i--) {
    z = fixpt[i];
    while (l < z) {
      l++;
      if (pno[l] > 0) {
        (*pno)++;
        pno[*pno] = pno[l];
      }
    }
    lorb[i] = orbitsv(base[i], svptr[i], 0);
  }
  fprintf(op, "%4d%4d%4d%4d\n", npt, *pno, nb, 3);
  printbaselo(nb, base, lorb);
  printpsv(nb, pno, svptr);
  for (i = 1; i <= *pno; i++)
    fprintf(op, "%4d", pno[i]);
  fprintf(op, "\n");
  return (0);
}
