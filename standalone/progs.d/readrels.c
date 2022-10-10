#include "defs.h"

#define MLEN 500
#define MDIM 51
#define MGENS 53
/* MLEN= max length of a relation.
   MDIM = max dimension of module, or no. of gens. of multiplier,
   MGENS=max no. of generators.
*/

short rel[MLEN], gch[128], gno[53], order[MDIM];
char  gap, inf1[80], inf2[80], io[80], outf[80], genletter[MGENS];
/* inf1=gpname.tc, inf2=gpname.inmat (if not -a or -m)
   io=gpname.subname then gpname.subname.rel. outf=gpname.subname.er if not -a
   By default subname=sg. Could be psg or nsg if P=G or N=G.
*/
FILE *ip, *ip2, *iop, *op;

static int readrel(int no);

static int inv(int n)
{
  if (n % 2 == 0)
    return (n + 1);
  return (n - 1);
}

static void inperr(int n)
{
  fprintf(stderr, "Input error in relation no %d\n", n);
}

static int digit(int j)
{
  if (j >= '0' && j <= '9')
    return (1);
  else
    return (0);
}

static int letter(int j)
{
  if ((j >= 'a' && j <= 'z') || (j >= 'A' && j <= 'Z'))
    return (1);
  else
    return (0);
}

static void snl_ip(void)
{
  while (getc(ip) != '\n')
    ;
}

static void snl_iop(void)
{
  while (getc(iop) != '\n')
    ;
}

int main(int argc, char * argv[])
{
  short i, j, k, l, n, np, nb, ng, nsg, nr, ch, dim, ngext, nrext, rno, ct;
  char  c, err, mult, append, split;
  int   arg;
  gap = err = append = mult = split = 0;
  arg = 1;
  while (argv[arg][0] == '-') {
    if (argv[arg][1] == 'a')
      append = 1;
    else if (argv[arg][1] == 's')
      split = 1;
    else if (argv[arg][1] == 'm')
      mult = 1;
    else if (argv[arg][1] == 'g')
      gap = 1;
    else {
      err = 1;
      goto error;
    }
    arg++;
  }
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcpy(inf1, argv[arg]);
  arg++;
  strcat(inf1, ".");
  strcpy(io, inf1);
  if (append == 0) {
    if (gap) {
      strcpy(outf, inf1);
      strcat(outf, "rvals");
    }
    else {
      strcpy(outf, inf1);
      strcat(outf, "ext.tc");
    }
    if (mult == 0) {
      strcpy(inf2, inf1);
      strcat(inf2, "inmat");
    }
  }
  strcat(inf1, "tc");
  if (argc <= arg)
    strcat(io, "sg");
  else
    strcat(io, argv[arg]);

  /* If appending, then first we read from gpname.sg to find nos of tc
     generators of G. These are on the last line of the file.
  */
  if (append) {
    if ((ip = fopen(io, "r")) == 0) {
      fprintf(stderr, "Cannot open file %s.\n", io);
      exit(1);
    }
    fscanf(ip, "%hd%hd%hd", &i, &np, &j);
    /* Go to last line of file */
    for (k = 1; k <= 3 + np + j; k++)
      snl_ip();
    for (i = 1; i <= 52; i++)
      gno[j] = -1;
    for (i = 1; i <= np; i++) {
      fscanf(ip, "%hd", &j);
      if (j <= 52)
        gno[j] = 2 * i - 2;
    }
    fclose(ip);
  }

  /* Now we open the Todd Coxeter file, and assign the above determined
     gen. nos. to the generators (which are letters).
  */
  if ((ip = fopen(inf1, "r")) == 0) {
    fprintf(stderr, "Cannot open file %s.\n", inf1);
    exit(1);
  }
  fscanf(ip, "%hd%hd%hd", &ng, &nsg, &nr);
  if (append)
    for (i = 1; i <= ng; i++) {
      if (gno[i] == -1) {
        fprintf(stderr, "Gen. no. %d does not occur in %s.\n", i, io);
        fprintf(stderr, "egrun should be run with -f option.\n");
        exit(1);
      }
    }
  for (i = 1; i < 128; i++)
    gch[i] = -1;
  for (i = 1; i <= ng; i++) {
    ch = ' ';
    while (ch == ' ' || ch == '\n')
      ch = getc(ip);
    if (letter(ch) == 0) {
      fprintf(stderr, "Generators must be letters.\n");
      exit(1);
    }
    if (gch[ch] != -1) {
      fprintf(stderr, "Repeated generator.\n");
      exit(1);
    }
    gch[ch] = gno[i];
    genletter[i] = ch;
  }
  snl_ip();
  for (i = 1; i <= nsg; i++)
    snl_ip();

  /* If appending, we open the file gpname.sg.rel for this purpose. Otherwise,
     we find the dimension of the normal subgroup, and write the obvious
     relations (using gpname.inmat if non-central) to the Todd-Coxeter type
     file.
  */

  if (append) {
    strcat(io, ".rel");
    if ((iop = fopen(io, "r")) == 0) {
      fprintf(stderr, "Cannot open file %s.\n", io);
      exit(1);
    }
    fscanf(iop, "%hd", &nb);
    fclose(iop);
    iop = fopen(io, "w");
    fprintf(iop, "%4d   %4d", nb, nr);
    for (i = 2; i <= nb; i++)
      fprintf(iop, "%4d", 0);
    fprintf(iop, "\n");
  }
  else {
    if (split == 0) {
      strcat(io, ".er");
      if ((iop = fopen(io, "r")) == 0) {
        fprintf(stderr, "Cannot open file %s.\n", io);
        exit(1);
      }
      fscanf(iop, "%hd%hd", &i, &dim);
      if (dim >= MDIM) {
        fprintf(stderr, "Too many generators of module. Increase MDIM.\n");
        exit(1);
      }
      snl_iop();
      fscanf(iop, "%hd", &k);
      snl_iop();
      if (mult) {
        for (i = 1; i <= dim; i++)
          fscanf(iop, "%hd", order + i);
        snl_iop();
      }
      if (k != nr) {
        fprintf(stderr, "No. of relations wrong in %s.\n", io);
        exit(1);
      }
    }
    if (mult == 0) {
      if ((ip2 = fopen(inf2, "r")) == 0) {
        fprintf(stderr, "Cannot open file %s.\n", inf2);
        exit(1);
      }
      fscanf(ip2, "%hd%hd%hd", &i, &j, &k);
      if (split) {
        dim = j;
        if (dim >= MDIM) {
          fprintf(stderr, "Too many generators of module. Increase MDIM.\n");
          exit(1);
        }
      }
      if ((j != dim) || (k != ng)) {
        fprintf(stderr, "Conflicting data in line 1 of %s.\n", inf2);
        exit(1);
      }
      for (k = 1; k <= dim; k++)
        order[k] = i;
    }

    op = fopen(outf, "w");
    if (gap)
      fprintf(op, "COHOMOLO.RelVals:=[\n");
    else {
      ngext = ng + dim;
      nrext = nr + dim * (dim + 1) / 2 + dim * ng;
      fprintf(op, "%4d%4d%4d\n", ngext, 0, nrext);
      for (i = 1; i <= ng; i++)
        putc(genletter[i], op);
      for (i = 1; i <= dim; i++)
        fprintf(op, "M[%d]", i);
      putc('\n', op);
      fprintf(stderr, "Relations of Extension are as follows:\n");
      fprintf(stderr, "     (These are being written to file %s.)\n", outf);
      for (i = 1; i <= dim; i++) {
        fprintf(op, "M[%d]%d\n", i, order[i]);
        fprintf(stderr, "M[%d]%d\n", i, order[i]);
      }
      for (j = 1; j < dim; j++)
        for (i = j + 1; i <= dim; i++) {
          fprintf(op, "M[%d]-M[%d]-M[%d]M[%d]\n", j, i, j, i);
          fprintf(stderr, "M[%d]-M[%d]-M[%d]M[%d]\n", j, i, j, i);
        }
      for (i = 1; i <= ng; i++)
        for (j = 1; j <= dim; j++) {
          fprintf(op, "%c-M[%d]%c", genletter[i], j, genletter[i]);
          fprintf(stderr, "%c-M[%d]%c", genletter[i], j, genletter[i]);
          if (mult) {
            fprintf(op, "M[%d]-\n", j);
            fprintf(stderr, "M[%d]-\n", j);
          }
          else {
            l = 0;
            for (k = 1; k <= dim; k++) {
              fscanf(ip2, "%hd", order + k);
              if (order[k] != 0)
                l++;
            }
            if (l > 1) {
              fprintf(op, "(");
              fprintf(stderr, "(");
            }
            for (k = 1; k <= dim; k++) {
              if ((n = order[k]) != 0) {
                if (n == 1) {
                  if (l == 1) {
                    fprintf(op, "M[%d]-\n", k);
                    fprintf(stderr, "M[%d]-\n", k);
                  }
                  else {
                    fprintf(op, "M[%d]", k);
                    fprintf(stderr, "M[%d]", k);
                  }
                }
                else {
                  if (l == 1) {
                    fprintf(op, "M[%d]%d\n", k, -n);
                    fprintf(stderr, "M[%d]%d\n", k, -n);
                  }
                  else {
                    fprintf(op, "M[%d]%d", k, n);
                    fprintf(stderr, "M[%d]%d", k, n);
                  }
                }
              }
            } /* k loop */
            if (l > 1) {
              fprintf(op, ")-\n");
              fprintf(stderr, ")-\n");
            }
          } /* mult==0 */
        }   /* i and j loop */
    }
    if (mult == 0)
      fclose(ip2);
  } /* append=0 */

  /* Now we are ready to read the relations, and write them on the end
     of the file gpname.sg.rel if appending, or write their values out
     to a Todd-Coxeter type file if constructing an extension.
  */
  for (rno = 1; rno <= nr; rno++) {
    if (append) {
      if (readrel(rno) == -1)
        exit(1);
      for (j = 0; j <= *rel; j++)
        fprintf(iop, "%4d", rel[j]);
      fprintf(iop, "\n");
    }
    else {
      if (gap == 0)
        while ((c = getc(ip)) != '\n') {
          putc(c, op);
          putc(c, stderr);
        }
      if (split == 0) {
        snl_iop();
        fscanf(iop, "%hd", &l);
        if (l > 2 && gap == 0) {
          putc('(', op);
          putc('(', stderr);
        }
        if (gap)
          putc('[', op);
        ct = 1;
        for (i = 1; i <= l / 2; i++) {
          fscanf(iop, "%hd%hd", &j, &k);
          if (gap) {
            if (j < ct) {
              fprintf(stderr, "Relation generators in wrong order.\n");
              exit(1);
            }
            while (ct < j) {
              putc('0', op);
              putc(',', op);
              ct++;
            }
            ct++;
            fprintf(op, "%d", k);
            if (j < dim)
              putc(',', op);
          }
          else {
            if (k == 1) {
              if (l == 2) {
                fprintf(op, "M[%d]-", j);
                fprintf(stderr, "M[%d]-", j);
              }
              else {
                fprintf(op, "M[%d]", j);
                fprintf(stderr, "M[%d]", j);
              }
            }
            else {
              if (l == 2) {
                fprintf(op, "M[%d]%d", j, -k);
                fprintf(stderr, "M[%d]%d", j, -k);
              }
              else {
                fprintf(op, "M[%d]%d", j, k);
                fprintf(stderr, "M[%d]%d", j, k);
              }
            }
          }
        }
        if (l > 2 && gap == 0) {
          fprintf(op, ")-");
          fprintf(stderr, ")-");
        }
        if (gap) {
          while (ct <= dim) {
            putc('0', op);
            if (ct < dim)
              putc(',', op);
            ct++;
          }
          putc(']', op);
          if (rno < nr)
            putc(',', op);
          putc('\n', op);
        }
        else {
          fprintf(op, "\n");
          fprintf(stderr, "\n");
        }
        snl_iop();
      } /* split==0 */
      else {
        fprintf(op, "\n");
        fprintf(stderr, "\n");
      }
    }
  }
  if (gap)
    fprintf(op, "];\n");
  fclose(ip);
  if (split == 0)
    fclose(iop);
  if (append == 0)
    fclose(op);

error:
  if (err) {
    fprintf(stderr, "Usage:  readrels [-g] [-a/-s/-m] gpname [sub].\n");
    exit(1);
  }
  exit(0);
}

int readrel(int no)
/* Reads relation number "no" into array rel, after
   expanding powers, brackets etc.
   Each relation is preceded by its length.
*/
{
  short stbr, endbr, exp, l, m, n, ptr, ch, len;
  char  gotg, br, clbr, emptybr;
  gotg = 0;
  br = 0;
  clbr = 0;
  emptybr = 0;
  stbr = 0;
  endbr = 0;
  ptr = 0;
  len = 0;
  ch = getc(ip);
  while (ch != '\n') {
    if (ch == ' ')
      ch = getc(ip);
    else if (ch == '(') {
      if (br || clbr) {
        inperr(no);
        return (-1);
      }
      emptybr = 1;
      br = 1;
      stbr = ptr + 1;
      gotg = 0;
      ch = getc(ip);
    }
    else if (ch == ')') {
      if (br == 0 || emptybr) {
        inperr(no);
        return (-1);
      }
      br = 0;
      gotg = 0;
      clbr = 1;
      endbr = ptr;
      ch = getc(ip);
    }
    else if (letter(ch)) {
      if (clbr || gch[ch] == -1) {
        inperr(no);
        return (-1);
      }
      emptybr = 0;
      gotg = 1;
      len++;
      ptr++;
      rel[ptr] = gch[ch];
      ch = getc(ip);
    }
    else if (ch == '-' || digit(ch)) {
      if (gotg == 0 && clbr == 0) {
        inperr(no);
        return (-1);
      }
      gotg = 0;
      if (ch == '-') {
        ch = getc(ip);
        if (digit(ch) == 0)
          exp = -1;
        else {
          exp = 0;
          while (digit(ch)) {
            exp *= 10;
            exp -= (ch - '0');
            ch = getc(ip);
          }
        }
      }
      else {
        exp = 0;
        while (digit(ch)) {
          exp *= 10;
          exp += (ch - '0');
          ch = getc(ip);
        }
      }
      if (exp == 0) {
        inperr(no);
        return (-1);
      }
      if (clbr) {
        if (exp < 0)
          for (m = stbr, n = endbr; m <= n; m++, n--) {
            if (m == n)
              rel[m] = inv(rel[m]);
            else {
              l = inv(rel[m]);
              rel[m] = inv(rel[n]);
              rel[n] = l;
            }
          }
        exp = abs(exp);
        exp--;
        clbr = 0;
        for (n = 1; n <= exp; n++)
          for (m = stbr; m <= endbr; m++) {
            len++;
            ptr++;
            rel[ptr] = rel[m];
          }
      }
      else {
        n = rel[ptr];
        if (exp < 0) {
          n++;
          rel[ptr] = n;
          exp = -exp;
        }
        exp--;
        for (m = 1; m <= exp; m++) {
          len++;
          ptr++;
          rel[ptr] = n;
        }
      }
    }
    else {
      inperr(no);
      return (-1);
    }
  }
  if (br || clbr) {
    inperr(no);
    return (-1);
  }
  rel[0] = len;
  return (0);
}
