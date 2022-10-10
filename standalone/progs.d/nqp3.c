#include "defs.h"

extern char  inf[], inf1[], outf1[];
extern short facexp, prime, exp, *rpf, *rpb, *eexpnt, **pcb, dim, onng,
    expnt[], **comptr[], *vec[], **mat[], cp[], *wf, *wc, **extno, **subno,
    chsdim, chpdim, gap;
extern FILE *ip, *op;

int comprels(void)
/* Used when crel is true to compute values of relators of P in the chosen
   stable extension of M by P.
*/
{
  short i, j, k, l, m, n, k1, *p, *q, *r, stabdim, *stabhom, **stabno, *orpf,
      nb, np, *rno, *covrel, **pgen, *ps, *pf, *v1, *v2, **cbmat;
  short sgn;
  stabdim = chpdim - chsdim;
  stabno = subno - stabdim;
  stabhom = rpf - 1;
  rpf += onng;
  orpf = rpf;
  for (i = 1; i <= stabdim; i++) {
    stabno[i] = rpf - 1;
    rpf += onng;
  }
  k = 0;
  for (i = 1; i <= onng; i++)
    if (subno[i] == 0 && extno[i] != 0) {
      k++;
      p = stabno[k];
      for (j = 1; j <= onng; j++)
        p[j] = 0;
      p[i] = 1;
      for (j = 1; j <= onng; j++)
        if ((q = subno[j]) != 0) {
          expand(q, rpf, onng);
          if ((l = rpf[i]) != 0)
            p[j] = prime - l;
        }
    }
  if (k != stabdim) {
    fprintf(stderr, "stabdim error.\n");
    return (-1);
  }
  if (gap == 0 && stabdim > 1)
    fprintf(stderr, "Basis of stable homomorphisms:\n");
  if (gap == 0 && stabdim > 1)
    for (i = 1; i <= stabdim; i++) {
      for (j = 1; j <= onng; j++)
        fprintf(stderr, "%3d", stabno[i][j]);
      fprintf(stderr, "\n");
    }
  for (i = 1; i <= onng; i++)
    stabhom[i] = 0;
  if (gap == 0 && stabdim > 1)
    fprintf(stderr,
            "Choose required stable hom as a vector in this basis!\n");
  if (gap)
  /* We just input the number of the basis of stable homs that we want */
  {
    scanf("%hd", &i);
    for (j = 1; j <= onng; j++)
      stabhom[j] = stabno[i][j];
  }
  else
    for (i = 1; i <= stabdim; i++) {
      if (stabdim > 1)
        scanf("%hd", rpf);
      else
        *rpf = 1;
      if (*rpf < 0) {
        for (j = 1; j <= onng; j++)
          scanf("%hd", stabhom + j);
        break;
      }
      if (*rpf != 0) {
        p = stabhom;
        r = p + onng;
        q = stabno[i] + 1;
        while (++p <= r) {
          *p += *q * (*rpf);
          *p %= prime;
          q++;
        }
      }
    }
  printf("Chosen hom is:\n");
  for (i = 1; i <= onng; i++)
    printf("%3d", stabhom[i]);
  printf("\n");
  rpf = orpf;
  strcpy(inf1, inf);
  strcat(inf1, "psgwds");
  ip = fopen(inf1, "r");
  if (ip == 0) {
    fprintf(stderr, "Cannot open file %s.\n", inf1);
    return (-1);
  }
  pgen = pcb;
  fscanf(ip, "%hd", &np);
  for (i = 1; i <= np; i++) {
    pgen[i] = rpf;
    fscanf(ip, "%hd", rpf);
    p = rpf;
    rpf += (1 + *p);
    while (++p < rpf)
      fscanf(ip, "%hd", p);
  }
  fclose(ip);
  /* Compute matrix to change basis of module back to the original, by
     using the base change matrices output by matcalc.
  */
  strcpy(inf1, inf);
  strcat(inf1, "cbmats");
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open file %s.\n", inf1);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd", &i, &j, &k);
  if (i != prime || j != dim || k != 2) {
    fprintf(stderr, "Error in line 1 of %s.\n", inf1);
    return (-1);
  }
  readmat(mat[1]);
  readmat(mat[2]);
  *cp = 2;
  cp[1] = 2;
  cp[2] = 1;
  prod(cp, mat[3]);
  inv(mat[3], mat[2]);
  trans(mat[2], mat[1]);
  cbmat = mat[1];
  v1 = mat[2][1];
  v2 = mat[3][1];
  fclose(ip);
  strcpy(inf1, inf);
  strcat(inf1, "psg.rel");
  ip = fopen(inf1, "r");
  if (ip == 0) {
    fprintf(stderr, "Cannot open file %s.\n", inf1);
    return (-1);
  }
  strcpy(outf1, inf);
  strcat(outf1, "psg.er");
  op = fopen(outf1, "w");
  fscanf(ip, "%hd", &nb);
  rno = rpf - 1;
  rpf += nb;
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%hd", rno + i);
  fprintf(op, "%4d%4d\n", nb, dim);
  for (i = 1; i <= nb; i++)
    fprintf(op, "%4d", rno[i]);
  fprintf(op, "\n");
  wf = rpf;
  for (i = 1; i <= rno[1]; i++) {
    fscanf(ip, "%hd", &l);
    covrel = rpb - l;
    p = covrel;
    while (++p <= rpb)
      fscanf(ip, "%hd", p);
    zero(expnt, eexpnt);
    for (j = l; j >= 1; j--) {
      wc = wf - 2;
      k = covrel[j];
      k1 = k / 2 + 1;
      m = *pgen[k1];
      if (k % 2 == 0) {
        sgn = 1;
        ps = pgen[k1] + 1;
        pf = ps + m - 2;
      }
      else {
        sgn = -1;
        pf = pgen[k1] + 1;
        ps = pf + m - 2;
      }
      while (1) {
        wc += 2;
        *wc = *ps;
        *(wc + 1) = *(ps + 1) * sgn;
        if (ps == pf)
          break;
        ps += (2 * sgn);
      }
      collect(wc, wf, 1);
    }
    fprintf(op, "%4d  ", l);
    for (j = 1; j <= l; j++)
      fprintf(op, "%4d", covrel[j]);
    fprintf(op, "\n");
    zero(v1, v1 + dim);
    for (n = 1; n <= exp; n++)
      if ((l = expnt[n]) != 0) {
        if (n <= facexp) {
          fprintf(stderr, "relation error. i,n,l=%d,%d,%d\n", i, n, l);
          return (-1);
        }
        for (j = 1; j <= dim; j++) {
          p = *(comptr[exp + j] + n);
          if (p != 0) {
            r = p + *p;
            while (++p < r)
              if ((k = stabhom[*p]) != 0) {
                v1[j] += (k * l * *(++p));
                v1[j] %= prime;
              }
              else
                ++p;
          }
        }
      }
    im(v1, v2, cbmat);
    l = 0;
    p = v2;
    while (++p <= v2 + dim)
      if (*p != 0)
        l += 2;
    fprintf(op, "%4d  ", l);
    p = v2;
    while (++p <= v2 + dim)
      if (*p != 0)
        fprintf(op, "%4d%4d", (int)(p - v2), *p);
    fprintf(op, "\n");
    if ((i == rno[1]) && (fscanf(ip, "%hd", &j) > 0)) {
      fprintf(op, "%4d\n", j);
      rno[1] += j;
    }
  } /* for (i=1;i<=rno[1];... */
  return (0);
}
