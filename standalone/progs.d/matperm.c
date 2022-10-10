#include "defs.h"

#define MSP 100000
#define MM 101
#define MV 1000
#define MPT 2000
#define MPR 2000
#define MDIM 300
#define MARG 3000

int   msp = MSP;
short prime, dim, *spv, **spm, nmat, ordim, npt, *fptr, *nfptr, **sp, **imsp,
    *endsp, mspace[MSP], *vec[MV], **mat[MM], pinv[MPR], *pptr[MM],
    *spptr[MPT + 1], mm = MM, mv = MV, mpt = MPT, mpr = MPR, mdim = MDIM,
                     marg = MARG;
char  orvec;
FILE *ip, *op;

int main(int argc, char * argv[])
{
  short nm, nv, i, j, pt, perm, *p1, *p2, *q1, arg;
  char  inf[80], outf[80], err;
  int   x;
  /* Defaults: inf=gpname.inmat  outf=gpname.inperm  */
  arg = 1;
  err = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcpy(inf, argv[1]);
  strcat(inf, ".");
  strcpy(outf, inf);
  arg++;
  if (argc <= arg)
    strcat(inf, "inmat");
  else
    strcat(inf, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf, "inperm");
  else
    strcat(outf, argv[arg]);
  if ((ip = fopen(inf, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf);
    exit(1);
  }
  fscanf(ip, "%hd%hd%hd", &prime, &dim, &nmat);
  if (prime > mpr) {
    fprintf(stderr, "prime too big. Increase MPR.\n");
    exit(1);
  }
  if (dim > mdim) {
    fprintf(stderr, "dim too big. Increase MDIM.\n");
    exit(1);
  }
  if (nmat > mm) {
    fprintf(stderr, "Too many mats. Increase MM.\n");
    exit(1);
  }
  setpinv();
  nm = nmat + 1;
  nv = (nm + 2) * dim;
  if (nv >= mv) {
    fprintf(stderr, "Too many vectors. Increase MV.\n");
    exit(1);
  }
  x = (nv + 1) * dim;
  if (x > msp - marg) {
    fprintf(stderr, "Out of space. Increase MSP.\n");
    exit(1);
  }
  /* Set matrix pointers as in mcp.c */
  for (i = 0; i <= nv; i++)
    vec[i] = mspace - 1 + i * dim;
  for (i = 0; i <= nmat; i++) {
    mat[i] = vec - 1 + i * dim;
    if (i > 0)
      readmat(mat[i]);
  }
  fclose(ip);
  endsp = mspace + msp - marg;
  spm = mat[0];
  spv = spm[1];
  sp = vec - 1 + nm * dim;
  imsp = sp + dim;
  fptr = vec[nv];
  for (i = 1; i <= nmat; i++) {
    pptr[i] = fptr;
    fptr += MPT;
  }
  if (fptr >= endsp) {
    fprintf(stderr, "Out of space. Increase MSP.\n");
    exit(1);
  }
reenter:
  printf("Input dim of orbit space (-1 for orbit of vector):    ");
  scanf("%hd", &ordim);
  if (ordim < -1 || ordim == 0 || ordim > dim) {
    printf("Inappropriate dimension.\n");
    while (getchar() != '\n')
      ;
    goto reenter;
  }
  else if (ordim == -1) {
    ordim = 1;
    orvec = 1;
  }
  else
    orvec = 0;
  printf("Input generating vectors of orbit space.\n");
  for (i = 1; i <= ordim; i++) {
    p1 = imsp[i];
    for (j = 1; j <= dim; j++)
      scanf("%hd", p1 + j);
  }
  npt = 1;
  if (orvec == 0)
    normalize();
  encode();
  spptr[1] = fptr;
  fptr = nfptr;
  if (fptr >= endsp) {
    fprintf(stderr, "Out of space. Increase MSP.\n");
    exit(1);
  }

  for (pt = 1; pt <= npt; pt++) {
    decode(pt);
    for (perm = 1; perm <= nmat; perm++) {
      for (j = 1; j <= ordim; j++)
        im(sp[j], imsp[j], mat[perm]);
      if (orvec == 0)
        normalize();
      encode();
      for (i = 1; i <= npt; i++) {
        p1 = fptr;
        p2 = nfptr;
        q1 = spptr[i] + 1;
        while (++p1 <= p2) {
          if (*p1 != *q1)
            break;
          q1++;
        }
        if (p1 > p2)
          break;
      }
      pptr[perm][pt] = i;
      if (i > npt) {
        npt++;
        if (npt > mpt) {
          fprintf(stderr, "Too many points.\n");
          exit(1);
        }
        spptr[npt] = fptr;
        fptr = nfptr;
        if (fptr >= endsp) {
          fprintf(stderr, "Out of space. Increase MSP.\n");
          exit(1);
        }
        if (npt % 10 == 0)
          printf("npt=%d.\n", npt);
      }
    }
  }
  op = fopen(outf, "w");
  fprintf(op, "%4d%4d%4d%4d\n", npt, nmat, 0, 0);
  for (i = 1; i <= nmat; i++)
    if (npt >= 1000) {
      for (j = 1; j <= npt; j++)
        fprintf(op, "%5d", pptr[i][j]);
      fprintf(op, "\n");
    }
    else {
      for (j = 1; j <= npt; j++)
        fprintf(op, "%4d", pptr[i][j]);
      fprintf(op, "\n");
    }
error:
  if (err) {
    fprintf(stderr, "Usage:    matperm gpname [inf] [outf]\n");
    exit(1);
  }
  exit(0);
}

int encode(void)
{
  short i, *p1, *p2, *q;
  char  c;
  nfptr = fptr;
  for (i = 1; i <= ordim; i++) {
    c = orvec;
    p1 = imsp[i];
    p2 = p1 + dim;
    q = p1;
    while (++q <= p2)
      if (*q != 0)
        if (c == 0) {
          c = 1;
          *(++nfptr) = q - p1;
        }
        else {
          *(++nfptr) = q - p1;
          *(++nfptr) = *q;
        }
    *(++nfptr) = 0;
  }
}

int decode(int n)
{
  short i, *p1, *p2, *ptr, *q;
  char  c;
  ptr = spptr[n];
  for (i = 1; i <= ordim; i++) {
    c = orvec;
    p1 = sp[i];
    p2 = p1 + dim;
    q = p1;
    while (++q <= p2)
      *q = 0;
    while (*(++ptr) != 0)
      if (c == 0) {
        p1[*ptr] = 1;
        c = 1;
      }
      else {
        p1[*ptr] = *(ptr + 1);
        ptr++;
      }
  }
}

int normalize(void)
{
  short rst, i, j, k, fac, *p1, *p2, *q1;
  int   n;
  rst = 1;
  for (i = 1; i <= dim; i++)
    for (j = rst; j <= ordim; j++)
      if (imsp[j][i] != 0) {
        p1 = imsp[j];
        if (j > rst) {
          imsp[j] = imsp[rst];
          imsp[rst] = p1;
        }
        fac = pinv[p1[i]];
        p2 = p1 + dim + 1;
        while (--p2 > p1) {
          n = *p2 * fac;
          n %= prime;
          *p2 = n;
        }
        for (k = 1; k <= ordim; k++) {
          if (imsp[k][i] == 0 || k == rst)
            continue;
          fac = prime - imsp[k][i];
          p2 = p1 + dim + 1;
          q1 = imsp[k] + dim;
          while (--p2 > p1) {
            n = *q1 + fac * *p2;
            n %= prime;
            *q1 = n;
            q1--;
          }
        }
        if (rst == ordim)
          return (0);
        rst++;
        break;
      }
}

int setpinv(void)
{
  int i, j;
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
