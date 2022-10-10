#include "defs.h"

extern char inf[], outf[], temp[], temp2[], mspace[], *vec[], **mat[], cvec[],
    pinv[], mstr[], full, intop, opt, aut;
extern short svec[], mdim, mv, mm, mpr, msl;
extern int   space;
char         prime;
short        dim, maxnull;

FILE *ip, *op;

int pkprog(void)
{
  char ngens, neg, fac, there, saved, gottm, gotsgn, err, digexp, nc, provn,
      trd, gotfac, heldc, len, **m0, **m1, *m01, irr;
  short nvecs, i, j, rno1, rno2, nlty, fnz, adno, sdim, qdim, tdim;
  int   msp;
  irr = opt ? 1 : 0;
  if ((ip = fopen(inf, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf);
    return (-1);
  }
  fscanf(ip, "%hd%hd%hd", &rno1, &dim, &rno2);
  prime = rno1;
  ngens = rno2;
  if (prime > mpr) {
    fprintf(stderr, "prime too big. Nothing to be done.\n");
    return (-1);
  }
  if (dim > mdim) {
    fprintf(stderr, "dim too big. Increase MDIM.\n");
    return (-1);
  }
  if (ngens >= mm) {
    fprintf(stderr, "Too many generators. Nothing to be done.\n");
    return (-1);
  }
  nvecs = dim * (ngens + 1);
  if (nvecs > mv) {
    fprintf(stderr, "Too many vectors. Increase MV.\n");
    return (-1);
  }
  msp = nvecs * dim;
  if (msp > space) {
    fprintf(stderr, "Not enough space. Increase SPACE.\n");
    return (-1);
  }
  seeknln();
  setpinv();
  /* Set vector and matrix pointers, as usual.
     We take care only to use ngens+1 matrices, thus saving as much storage as
     possible. This means using temporary files to store matrices
     occasionally.
  */
  for (i = 0; i < nvecs; i++)
    vec[i] = mspace - 1 + i * dim;
  for (i = 0; i <= ngens; i++)
    mat[i] = vec - 1 + i * dim;
  for (i = 1; i <= ngens; i++)
    readmat(mat[i]);
  fclose(ip);
  if (aut) {
    printf("Enter maximum value of nullity that should be handled:   ");
    scanf("%hd", &maxnull);
  }
  m0 = mat[0];
  m1 = mat[1];
  m01 = m0[1];
  /* The next section reads in the given expression in the generators, and
     computes the corresponding matrix as m0. A temporary file is used to
     remember the sum so far, whilst a product is being computed.
  */
  while (1) {
    printf(
        "Enter element of group algebra. End expression with $. Example:\n");
    printf("Enter 3a[3]a[3]a[2] + a[2] - a[4]  as  3.3*3*2 + 2 - 4$\n");
    fflush(stdout);
    there = 0;
    saved = 0;
    heldc = 0;
    err = 0;
    while (heldc != '$') {
      gotfac = 0;
      gotsgn = 0;
      gottm = 0;
      provn = 0;
      digexp = 0;
      len = 0;
      while (gottm == 0) {
        if (heldc) {
          nc = heldc;
          heldc = 0;
        }
        else
          nc = findc();
        if (gotsgn == 0) {
          if (nc == '-') {
            neg = 1;
            gotsgn = 1;
            digexp = 1;
            continue;
          }
          else if (nc == '+') {
            neg = 0;
            gotsgn = 1;
            digexp = 1;
            continue;
          }
          else if (there == 0) {
            neg = 0;
            gotsgn = 1;
            digexp = 1;
          }
        }
        if ((nc == '+' || nc == '-' || nc == '$') && digexp == 0 &&
            (len != 0 || provn != 0)) {
          heldc = nc;
          gottm = 1;
          if (gotfac == 0) {
            if (provn) {
              fac = neg ? prime - 1 : 1;
              len = 1;
              mstr[1] = provn - '0';
            }
          }
        }
        else if (nc == '.' && digexp == 0 && provn != 0 && gotfac == 0) {
          gotfac = 1;
          fac = provn - '0';
          fac %= prime;
          provn = 0;
          if (neg)
            fac = prime - fac;
          digexp = 1;
        }
        else if (nc == '*' && digexp == 0 && (provn != 0 || len != 0)) {
          if (gotfac == 0) {
            gotfac = 1;
            fac = neg ? prime - 1 : 1;
            len = 1;
            mstr[1] = provn - '0';
          }
          digexp = 1;
        }
        else if (digit(nc) && digexp) {
          if (gotfac == 0)
            provn = nc;
          else {
            len++;
            mstr[len] = nc - '0';
          }
          digexp = 0;
        }
        else {
          fprintf(stderr, "Format error; Try again.\n");
          err = 1;
          while (getchar() != '\n')
            ;
          break;
        }
      }
      if (err)
        break;
      if (len == 1) {
        if (there)
          sum(fac, mat[mstr[1]], m0);
        else
          ncopy(fac, mat[mstr[1]], m0);
      }
      else {
        if (there)
        /* We have a partial sum, and now another product to compute, so we
           must store the partial sum in a temporary file.
        */
        {
          op = fopen(temp, "w");
          printmat(m0);
          saved = 1;
          fclose(op);
        }
        ncopy(fac, mat[mstr[1]], m0);
        for (i = 2; i <= len; i++)
          prod(m0, mat[mstr[i]], m0);
        if (saved) {
          ip = fopen(temp, "r");
          for (i = 1; i <= dim; i++)
            rvecsum(1, m0[i]);
          fclose(ip);
          saved = 0;
        }
      }
      there = 1;
    }
    if (err)
      continue;
    /* err means the expression input was garbled */
    if (intop)
      for (i = 1; i <= dim; i++) {
        for (j = 1; j <= dim; j++)
          printf("%3d", m0[i][j]);
        printf("\n");
      }
    /* Store the matrix m0 in a temporary file temp. Then compute its nullity.
       (nlty writes generators of the null space to another temporary file,
       temp)
    */
    op = fopen(temp2, "w");
    printmat(m0);
    fclose(op);
    op = fopen(temp, "w");
    nlty = null(m0, m1);
    fclose(op);
    ip = fopen(inf, "r");
    seeknln();
    readmat(m1);
    fclose(ip);
    if (nlty == 0 || (aut && nlty > maxnull))
      continue;
    if (aut)
      break;
    printf("Shall we use this matrix (y/n)?  ");
    nc = findc();
    if (nc == 'y')
      break;
  }
  /* Now we have chosen the matrix that we are going to use for the
     irreducibility test. Next we arrange to consider the one-dimensional
     subspaces of the null-space in order. cvec is used for this purpose.
  */
  for (i = 1; i <= nlty; i++)
    cvec[i] = 0;
  fnz = nlty;
  adno = nlty;
  trd = 0;
  while (fnz != 0) {
    cvec[adno]++;
    for (i = 1; i <= dim; i++)
      m01[i] = 0;
    /* Compute the current vector in the nullspace by reading from the file
       temp, and store it in m01
    */
    ip = fopen(temp, "r");
    for (i = 1; i <= nlty; i++) {
      if ((fac = cvec[i]) == 0)
        seeknln();
      else {
        rvecsum(fac, m01);
        seeknln();
      }
    }
    fclose(ip);
    if (intop) {
      for (j = 1; j <= dim; j++)
        printf("%3d", m01[j]);
      printf("\n");
      if (opt) {
        printf("Shall we skip this one?   ");
        nc = findc();
        if (nc == 'y')
          goto next;
      }
    }
    /* Now compute the subspace it generates */
    sdim = spgen(m0, ngens);
    if (intop) {
      for (i = 1; i <= dim; i++)
        printf("%3d", svec[i]);
      printf("\n");
      for (i = 1; i <= dim; i++) {
        for (j = 1; j <= dim; j++)
          printf("%3d", m0[i][j]);
        printf("\n");
      }
    }
    if (sdim < dim)
      /* sdim<dim means the module is reducible */
      if (opt == 0)
        break;
      else {
        printf("Shall we compute (y/n)?  ");
        nc = findc();
        if (nc == 'y') {
          opt = 0;
          break;
        }
        irr = 0;
      }
  next:
    /* The last vector generated the whole space, so go on to the next */
    adno = nlty;
    while (cvec[adno] == prime - 1) {
      cvec[adno] = 0;
      adno--;
    }
    if (adno <= fnz) {
      if (adno == fnz) {
        cvec[adno] = 0;
        adno--;
      }
      fnz--;
    }
  }
  if (opt || sdim == dim)
  /* Each vector in the null-space generated the whole space. Finally we have
     to transpose everything, and try just one vector in the null-space.
  */
  {
    printf("Transforming.\n");
    fflush(stdout);
    trd = 1;
    ip = fopen(temp2, "r");
    readmat(m0);
    fclose(ip);
    trans(m0, m1);
    op = fopen(temp, "w");
    null1(m1, m0);
    fclose(op);
    /* That ruined m1 (the first generator), so we must re-read it */
    ip = fopen(inf, "r");
    seeknln();
    readmat(m1);
    fclose(ip);
    for (i = 1; i <= ngens; i++) {
      trans(mat[i], m0);
      copy(m0, mat[i]);
    }
    ip = fopen(temp, "r");
    for (i = 1; i <= dim; i++) {
      fscanf(ip, "%hd", &rno1);
      m01[i] = rno1;
    }
    fclose(ip);
    sdim = spgen(m0, ngens);
    if (opt && sdim < dim) {
      printf("Shall we compute (y/n)?  ");
      nc = findc();
      if (nc == 'y')
        opt = 0;
      else
        irr = 0;
    }
    if (intop) {
      for (i = 1; i <= dim; i++)
        printf("%3d", svec[i]);
      printf("\n");
      for (i = 1; i <= dim; i++) {
        for (j = 1; j <= dim; j++)
          printf("%3d", m0[i][j]);
        printf("\n");
      }
    }
  }

  if (opt && irr == 0) {
    unlink(temp);
    unlink(temp2);
    return (0);
  }
  if (sdim == dim) {
    printf("Module is irreducible.\n");
    unlink(temp);
    unlink(temp2);
    return (0);
  }
  /* Space is reducible. Finally, we compute the generators of the subspace
     and quotient space
  */
  qdim = dim - sdim;
  strcpy(outf, inf);
  if (full) {
    strcat(outf, "f");
    tdim = dim;
  }
  else {
    strcat(outf, "s");
    tdim = sdim;
  }
  /* if we have transposed, then we first output to the temporary file,
     and then later read it back in and transpose it back again.
  */
  if (trd)
    op = fopen(temp, "w");
  else {
    op = fopen(outf, "w");
    fprintf(op, "%3d%5d%3d\n", prime, tdim, ngens);
  }
  opnmat(m0, ngens, tdim, 1);
  if (trd == 0 || full)
    fclose(op);
  if (full == 0) {
    strcpy(outf, inf);
    strcat(outf, "q");
    if (trd == 0) {
      op = fopen(outf, "w");
      fprintf(op, "%3d%5d%3d\n", prime, qdim, ngens);
    }
    opnmat(m0, ngens, dim, sdim + 1);
    fclose(op);
  }
  if (trd) {
    if (full == 0)
      dim = sdim;
    /* Retranspose */
    ip = fopen(temp, "r");
    op = fopen(outf, "w");
    fprintf(op, "%3d%5d%3d\n", prime, dim, ngens);
    for (i = 1; i <= ngens; i++) {
      readmat(m0);
      trans(m0, m1);
      printmat(m1);
    }
    fclose(op);
    if (full == 0) {
      dim = qdim;
      strcpy(outf, inf);
      strcat(outf, "s");
      op = fopen(outf, "w");
      fprintf(op, "%3d%5d%3d\n", prime, dim, ngens);
      for (i = 1; i <= ngens; i++) {
        readmat(m0);
        trans(m0, m1);
        printmat(m1);
      }
      fclose(ip);
      fclose(op);
    }
  }
  unlink(temp);
  unlink(temp2);
  return (0);
}

void seeknln(void)
{
  while (getc(ip) != '\n')
    ;
}

void setpinv(void)
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

int findc(void)
{
  char c;
  while ((c = getchar()) == ' ' || c == '\n')
    ;
  return (c);
}

int digit(int c)
{
  if (c >= '1' && c <= '9')
    return (1);
  else
    return (0);
}
