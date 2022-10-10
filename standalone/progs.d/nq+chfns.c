#include "defs.h"

#define MEXP 101
#define RSP 5000
#define MCL 21
#define PTRSP 1000
#define MM 50

/*  exponent, RSP=basic space, PTRSP=basic pointer space, MCL=max class, */

char inf3[80], inf4[80], outf[80];
/* Always: inf3=gpname.pcp (Output of pcrun)
           inf4=gpname.pgmat
           outf=mats
*/
extern char inf1[];
short mexp = MEXP, mcl = MCL, no, rel[RSP], wt[MEXP], exp, *rpf, *rpb, **pcb,
      d1[MEXP], d2[MEXP], *pcptr[PTRSP], **powptr[MEXP], **comptr[MEXP],
      rsp = RSP, ptrsp = PTRSP, maxv, maxm, ngens;
extern short prime, dim, *spv, **spm, mspace[], *vec[], **mat[], cp[], pinv[],
    opmats, mm, mv;
extern int msp;

int calcmats(void)
{
  short i, j, k, ct, *p, *q, **swop;
  strcpy(inf3, inf1);
  strcat(inf3, "pcp");
  strcpy(inf4, inf1);
  strcat(inf4, "pgmat");
  if (opmats) {
    strcpy(outf, inf1);
    strcat(outf, "mats");
  }
  if (rdmats() == -1)
    return (-1);
  if (ingp() == -1)
    return (-1);
  ct = ngens;
  printf("Computing matrices.\n");
  /* Matrices for all pcp gens of P  are now computed */
  if (maxm < 3 * exp + 2) {
    fprintf(stderr, "Not enough mat space. Increase MSP (of MV or MM).\n");
    return (-1);
  }
  for (i = exp; i >= 1; i--)
    if (wt[i] == 1) {
      if (i > ct) {
        swop = mat[i];
        mat[i] = mat[ct];
        mat[ct] = swop;
        for (j = 1; j <= dim; j++)
          if (d2[exp + j] == ct)
            d2[exp + j] = i;
      }
      inv(mat[i], mat[exp + i]);
      ct--;
    }
  if (ct != 0) {
    fprintf(stderr, "No of pgens wrong.\n");
    return (-1);
  }
  for (i = 2; i <= exp; i++)
    if (wt[i] > 1) {
      p = (d1[i] == d2[i]) ? *powptr[d1[i]] : *(comptr[d1[i]] + d2[i]);
      q = p + *p - 2;
      *cp = 0;
      while (--q > p) {
        k = *(q + 1);
        for (j = 1; j <= k; j++)
          cp[++(*cp)] = *q + exp;
        q--;
      }
      if (d1[i] == d2[i])
        for (j = 1; j <= prime; j++)
          cp[++(*cp)] = d1[i];
      else {
        cp[++(*cp)] = d1[i] + exp;
        cp[++(*cp)] = d2[i] + exp;
        cp[++(*cp)] = d1[i];
        cp[++(*cp)] = d2[i];
      }
      prod(cp, mat[i]);
      inv(mat[i], mat[i + exp]);
    }
  for (i = 1; i <= exp; i++)
    trans(mat[i + exp], mat[i]);
  if (opmats) {
    FILE * op = fopen(outf, "w");
    fprintf(op, "%4d%4d%4d\n", prime, dim, exp);
    for (i = 1; i <= exp; i++)
      printmat(mat[i]);
    fclose(op);
  }
  return (0);
}

int rdmats(void)
/* reads matrices of generators of P */
{
  short  i;
  int    quot;
  FILE * ip;
  ip = fopen(inf4, "r");
  if (ip == 0) {
    fprintf(stderr, "Cannot open %s\n ", inf4);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd", &prime, &dim, &ngens);
  setpinv();
  quot = msp / dim;
  if (quot > mv)
    quot = mv;
  maxv = quot;
  for (i = 0; i < maxv; i++)
    vec[i] = mspace - 1 + i * dim;
  maxm = maxv / dim;
  if (maxm >= mm)
    maxm = mm - 1;
  for (i = 0; i <= maxm; i++)
    mat[i] = vec - 1 + i * dim;
  spm = mat[0];
  spv = spm[1];
  if (maxm < ngens) {
    fprintf(stderr, "Not enough mat space. Increase MSP (of MV or MM).\n");
    return (-1);
  }
  for (i = 1; i <= 4; i++)
    while (getc(ip) != '\n')
      ;
  for (i = 1; i <= ngens; i++)
    readmat(mat[i]);
  fclose(ip);
  return (0);
}

int ingp(int inp)
/* Read in output of respcrun -s */
{
  short  i, j, k, l, m, *orpf, **pcp;
  FILE * ip;
  ip = fopen(inf3, "r");
  if (ip == 0) {
    fprintf(stderr, "Cannot open %s\n", inf3);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd%hd%hd%hd", &prime, &exp, &i, &no, &j, &m);
  if (exp >= mexp) {
    fprintf(stderr, "exp too big. Increase MEXP.\n");
    return (-1);
  }
  if (m != 0) {
    fprintf(stderr, "Wrong version of nq.\n");
    return (-1);
  }
  for (i = 1; i <= exp; i++)
    fscanf(ip, "%hd", wt + i);
  for (i = 1; i <= exp; i++)
    fscanf(ip, "%hd", d1 + i);
  for (i = 1; i <= exp; i++)
    fscanf(ip, "%hd", d2 + i);
  rpf = rel;
  rpb = rel + rsp - 1;
  pcp = pcptr;
  for (i = 2; i <= no; i++) {
    comptr[i] = pcp - 1;
    for (j = 1; j < i; j++) {
      fscanf(ip, "%hd%hd", &k, &l);
      if (l == 0)
        *(pcp++) = 0;
      else {
        *(pcp++) = rpf + 1;
        *(rpf++) = k;
        *rpf = l;
        orpf = rpf;
        rpf += (l + 1);
        while ((++orpf) < rpf)
          fscanf(ip, "%hd", orpf);
      }
    }
  }
  for (i = 1; i <= no; i++) {
    powptr[i] = pcp;
    fscanf(ip, "%hd%hd", &k, &l);
    if (l == 0)
      *(pcp++) = 0;
    else {
      *(pcp++) = rpf + 1;
      *(rpf++) = k;
      *rpf = l;
      orpf = rpf;
      rpf += (l + 1);
      while ((++orpf) < rpf)
        fscanf(ip, "%hd", orpf);
    }
  }
  for (i = no + 1; i <= exp; i++) {
    comptr[i] = pcp - 1;
    for (j = 1; j < i; j++)
      *(pcp++) = 0;
    powptr[i] = pcp - 1;
    *(pcp++) = 0;
  }
  if (pcp - pcptr > ptrsp) {
    fprintf(stderr, "Not enough pointer space. Increase PTRSP.\n");
    return (-1);
  }
  if (rpf - rel > rsp) {
    fprintf(stderr, "Not enough space. Increase RSP.\n");
    return (-1);
  }
  fclose(ip);
  return (0);
}

int setpinv(void)
{
  short i, j;
  int   sum;
  for (i = 0; i < prime; i++)
    pinv[i] = 0;
  for (i = 1; i < prime; i++)
    if (pinv[i] == 0)
      for (j = 1; j < prime; j++)
        if ((sum = i * j) % prime == 1) {
          pinv[i] = j;
          pinv[j] = i;
          break;
        }
}
