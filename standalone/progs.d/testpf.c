#include "defs.h"

#define PSP 60000
#define MP 200
#define NPT 8000
#define MB 50
char  opt[40], sysstring[40];
short mp = MP, mpt = NPT, mb = MB, npt, perm[PSP], cp[100], *pptr[MP],
      orb[NPT + 1], orno[NPT + 1], pno[MP], base[MB];
int   psp = PSP;
FILE *ip, *op;

int main(void)
{
  short i, j, m, n, y, np, nb, mxp, ord, im, len, *p, *q;
  int   quot;
  int   chct;
  char  nt, f, id, eq, fault, flnm[80];
reenter:
  printf("Do you wish to read permutations from a file (y/n)?  ");
  if (getchar() == 'y') {
    f = 1;
    snl();
    printf("Input filename:    ");
    scanf("%s", flnm);
    snl();
    if ((ip = fopen(flnm, "r")) == 0) {
      fprintf(stderr, "Cannot open %s.\n", flnm);
      goto reenter;
    }
    fscanf(ip, "%hd%hd%hd%hd", &npt, &np, &nb, &m);
    seeknln();
    if (nb != 0) {
      if (nb >= mb) {
        fprintf(stderr, "Too many base points. Increase MB.\n");
        return (-1);
      }
      for (i = 1; i <= nb; i++)
        fscanf(ip, "%hd", base + i);
      seeknln();
    }
    if (m != 0)
      seeknln();
  }
  else {
    f = 0;
    snl();
    printf("Input npt:   ");
    scanf("%hd", &npt);
    snl();
    nb = 0;
  }
  if (npt > mpt) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    exit(1);
  }
  quot = psp / npt;
  if (quot > mp)
    mxp = mp;
  else
    mxp = quot;
  printf("Perm nos 1 - %d can be read. Perm no 0 is always the identity\n",
         mxp - 1);
  for (i = 0; i < mxp; i++) {
    pptr[i] = perm - 1 + npt * i;
    pno[i] = 0;
  }
  p = pptr[0];
  for (i = 1; i <= npt; i++)
    p[i] = i;
  pno[0] = 1;
  if (f) {
    if (np >= mxp) {
      fprintf(stderr, "Not enough perm space to read perms from file.");
      fprintf(stderr, " Increase PSP (or MP).\n");
      exit(1);
    }
    for (i = 1; i <= np; i++) {
      if (readperm(pptr[i]) == 2) {
        fprintf(stderr, "Perm no %i is not a permutation.\n");
        exit(1);
      }
      pno[i] = 1;
      seeknln();
    }
    fclose(ip);
  }

  while (1) {
    printf("Choose option. Print 'l' for list.\n");
    scanf("%s", opt);
    if (strcmp(opt, "l") == 0) {
      snl();
      printf("rp  n                = read perm no n.\n");
      printf("pp  n                = print perm no n.\n");
      printf("ppc  n               = print perm no n in cycles.\n");
      printf("inv m n              = calc inverse p[n] of p[m].\n");
      printf("pr n l i(1)...i(l)   = calc prod p[n] of p[i(1)]...p[i(l)].\n");
      printf("im m x               = print image of point x under p[m].\n");
      printf("fp n                 = print fixed pts of p[n].\n");
      printf("fpn l i(1)...i(l)    = print common fixed pts of "
             "p[i(1)],...,p[i(l)].\n");
      printf("cyc n x              = print cycle of point x under p[n].\n");
      printf("ord n                = print order of p[n].\n");
      printf("orb l i(1)...i(l)    = print orbits of p[i(1)]...p[i(l)].\n");
      printf("rb                   = read in a base for the group.\n");
      printf("op filename          = output some perms to 'filename'.\n");
      printf("opord filename       = output some perms to 'filename' and\n");
      printf(
          "                       compute order of group they generate.\n");
      printf("                      (Only works if base known for G!).\n");
      printf(
          "rf n filename        = read in additional perms from filename\n");
      printf("                       starting at perm. no. n.\n");
      printf("q                    = quit.\n");
    }
    else if (strcmp(opt, "rp") == 0) {
      scanf("%hd", &n);
      snl();
      if (n <= 0 || n >= mxp)
        fprintf(stderr, "Invalid perm.no %d.\n", n);
      else {
        while (rp(pptr[n]) == 2) {
          fprintf(stderr, "That is not a permutation. Try again!\n");
          snl();
        }
        pno[n] = 1;
      }
    }
    else if (strcmp(opt, "pp") == 0) {
      scanf("%hd", &n);
      snl();
      if (n < 0 || n >= mxp || pno[n] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", n);
        continue;
      }
      p = pptr[n];
      for (i = 1; i <= npt; i++)
        printf("%4d", p[i]);
      printf("\n");
    }
    else if (strcmp(opt, "ppc") == 0) {
      scanf("%hd", &n);
      snl();
      if (n < 0 || n >= mxp || pno[n] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", n);
        continue;
      }
      for (i = 1; i <= npt; i++)
        orb[i] = 1;
      p = pptr[n];
      id = 1;
      chct = 0;
      for (m = 1; m <= npt; m++)
        if (orb[m]) {
          if ((im = p[m]) != m) {
            id = 0;
            orb[im] = 0;
            printf("(%d,%d", m, im);
            chct += (num_digits(m) + num_digits(im) + 2);
            while ((im = p[im]) != m) {
              if (chct >= 75) {
                printf(",\n%d", im);
                chct = num_digits(im);
              }
              else {
                printf(",%d", im);
                chct += (1 + num_digits(im));
              }
              orb[im] = 0;
            }
            printf(")");
            chct++;
            if (chct >= 72) {
              printf("\n");
              chct = 0;
            }
          }
        }
      if (id)
        printf("Perm no %d is the identity.\n", n);
      else
        printf("\n");
    }
    else if (strcmp(opt, "inv") == 0) {
      scanf("%hd%hd", &m, &n);
      snl();
      if (m < 0 || m >= mxp || pno[m] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", m);
        continue;
      }
      if (n <= 0 || n >= mxp) {
        fprintf(stderr, "Invalid perm. no. %d.\n", n);
        continue;
      }
      invert(pptr[m], pptr[n]);
      pno[n] = 1;
    }
    else if (strcmp(opt, "pr") == 0) {
      scanf("%hd%hd", &m, cp);
      if (m <= 0 || m >= mxp) {
        fprintf(stderr, "Invalid perm no. %d\n", m);
        snl();
        continue;
      }
      fault = 0;
      for (i = 1; i <= *cp; i++) {
        scanf("%hd", cp + i);
        n = *(cp + i);
        if (n < 0 || n >= mxp || pno[n] == 0) {
          fprintf(stderr, "Perm no %d is not defined.\n", n);
          fault = 1;
          snl();
          break;
        }
      }
      if (fault)
        continue;
      snl();
      p = pptr[m];
      for (i = 1; i <= npt; i++)
        p[i] = image(i);
      pno[m] = 1;
    }
    else if (strcmp(opt, "im") == 0) {
      scanf("%hd%hd", &n, &i);
      snl();
      if (n < 0 || n >= mxp || pno[n] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", n);
        continue;
      }
      if (i <= 0 || i > npt) {
        fprintf(stderr, "Inappropriate point %d.\n", i);
        continue;
      }
      printf("%d\n", pptr[n][i]);
    }
    else if (strcmp(opt, "fp") == 0) {
      scanf("%hd", &n);
      snl();
      if (n < 0 || n >= mxp || pno[n] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", n);
        continue;
      }
      p = pptr[n];
      nt = 1;
      for (i = 1; i <= npt; i++)
        if (p[i] == i) {
          nt = 0;
          printf("%4d", i);
        }
      if (nt)
        printf("No fixed points.\n");
      else
        printf("\n");
    }
    else if (strcmp(opt, "fpn") == 0) {
      scanf("%hd", cp);
      fault = 0;
      for (i = 1; i <= *cp; i++) {
        scanf("%hd", cp + i);
        n = *(cp + i);
        if (n < 0 || n >= mxp || pno[n] == 0) {
          fprintf(stderr, "Perm no %d is not defined.\n", n);
          fault = 1;
          snl();
          break;
        }
      }
      if (fault)
        continue;
      snl();
      nt = 1;
      for (n = 1; n <= npt; n++) {
        f = 1;
        for (i = 1; i <= *cp; i++)
          if (pptr[cp[i]][n] != n) {
            f = 0;
            break;
          }
        if (f) {
          nt = 0;
          printf("%4d", n);
        }
      }
      if (nt)
        printf("No fixed points.\n");
      else
        printf("\n");
    }
    else if (strcmp(opt, "cyc") == 0) {
      scanf("%hd%hd", &n, &i);
      snl();
      if (i <= 0 || i > npt) {
        fprintf(stderr, "Inappropriate point %d.\n", i);
        continue;
      }
      if (n < 0 || n >= mxp || pno[n] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", n);
        continue;
      }
      p = pptr[n];
      if ((im = p[i]) == i)
        printf("%d is fixed under p[%d].\n", i, n);
      else {
        printf("(%d,%d", i, im);
        while ((im = p[im]) != i)
          printf(",%d", im);
        printf(")\n");
      }
    }
    else if (strcmp(opt, "ord") == 0) {
      scanf("%hd", &n);
      snl();
      if (n < 0 || n >= mxp || pno[n] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", n);
        continue;
      }
      p = pptr[n];
      id = 1;
      for (i = 1; i <= npt; i++)
        orb[i] = 1;
      ord = 1;
      for (m = 1; m <= npt; m++)
        if (orb[m]) {
          if ((im = p[m]) != m) {
            id = 0;
            orb[im] = 0;
            len = 2;
            while ((im = p[im]) != m) {
              orb[im] = 0;
              len++;
            }
            ord = lcm(ord, len);
          }
        }
      if (id)
        printf("Perm no %d is the identity.\n", n);
      else
        printf("%d\n", ord);
    }
    else if (strcmp(opt, "orb") == 0) {
      scanf("%hd", cp);
      fault = 0;
      for (i = 1; i <= *cp; i++) {
        scanf("%hd", cp + i);
        n = *(cp + i);
        if (n < 0 || n >= mxp || pno[n] == 0) {
          fprintf(stderr, "Perm no %d is not defined.\n", n);
          fault = 1;
          snl();
          break;
        }
      }
      if (fault)
        continue;
      snl();
      allorbs();
    }
    else if (strcmp(opt, "eq") == 0) {
      scanf("%hd%hd", &m, &n);
      snl();
      if (n < 0 || n >= mxp || pno[n] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", n);
        continue;
      }
      if (m < 0 || m >= mxp || pno[m] == 0) {
        fprintf(stderr, "Perm no %d is not defined.\n", m);
        continue;
      }
      p = pptr[m];
      q = pptr[n];
      eq = 1;
      for (i = 1; i <= npt; i++)
        if (p[i] != q[i]) {
          eq = 0;
          break;
        }
      if (eq)
        printf("Permutations are equal.\n");
      else
        printf("Permutations are not equal on point %d.\n", i);
    }
    else if (strcmp(opt, "rb") == 0) {
      snl();
      printf("Input base points, preceded by their number.\n");
      scanf("%hd", &nb);
      if (nb >= mb) {
        fprintf(stderr, "Too many base points. Increase MB.\n");
        return (-1);
      }
      for (i = 1; i <= nb; i++)
        scanf("%hd", base + i);
      snl();
    }
    else if ((y = strcmp(opt, "opord")) == 0 || strcmp(opt, "op") == 0) {
      if (y == 0 && nb == 0) {
        fprintf(stderr, "Order can only be computed when base is known.\n");
        continue;
      }
      scanf("%s", flnm);
      snl();
      op = fopen(flnm, "w");
      printf("Input perm nos to be saved, preceded by their number.\n");
      scanf("%hd", &n);
      fprintf(op, "%3d %4d%4d%4d\n", npt, n, nb, 0);
      if (nb != 0)
        if (npt >= 10000) {
          for (i = 1; i <= nb; i++)
            fprintf(op, "%6d", base[i]);
          fprintf(op, "\n");
        }
        else if (npt >= 1000) {
          for (i = 1; i <= nb; i++)
            fprintf(op, "%5d", base[i]);
          fprintf(op, "\n");
        }
        else {
          for (i = 1; i <= nb; i++)
            fprintf(op, "%4d", base[i]);
          fprintf(op, "\n");
        }
      f = 0;
      for (i = 1; i <= n; i++) {
        scanf("%hd", &m);
        if (m < 0 || m >= mxp || pno[m] == 0) {
          printf("Perm no %d is not defined.\n", m);
          fclose(op);
          snl();
          unlink(flnm);
          f = 1;
          break;
        }
        printvec(pptr[m], 0);
      }
      if (f)
        continue;
      snl();
      fclose(op);
      if (y == 0) {
        sprintf(sysstring, "cp %s optxxx.ip", flnm);
        system(sysstring);
        system("gprun -b optxxx ip ip");
        unlink("optxxx.ip");
      }
    }
    else if (strcmp(opt, "rf") == 0) {
      scanf("%hd", &n);
      scanf("%s", flnm);
      snl();
      if ((ip = fopen(flnm, "r")) == 0) {
        fprintf(stderr, "Cannot open %s.\n", flnm);
        continue;
      }
      fscanf(ip, "%hd%hd%hd%hd", &i, &np, &j, &m);
      seeknln();
      if (i != npt) {
        fprintf(stderr, "npt does not agree.\n");
        fclose(ip);
        continue;
      }
      if (j != 0)
        seeknln();
      if (m != 0)
        seeknln();
      if (n <= 0 || np + n - 1 >= mxp) {
        fprintf(stderr, "Not enough perm space to read perms from file.");
        fclose(ip);
        continue;
      }
      for (i = 1; i <= np; i++) {
        if (readperm(pptr[i + n - 1]) == 2) {
          fprintf(stderr, "Perm no %i is not a permutation.\n");
          pno[i + n - 1] = 0;
          break;
        }
        pno[i + n - 1] = 1;
        seeknln();
      }
      fclose(ip);
    }
    else if (strcmp(opt, "q") == 0) {
      snl();
      break;
    }
    else {
      printf("Invalid option.\n");
      snl();
    }
  }
  exit(0);
}

int snl(void)
{
  while (getchar() != '\n')
    ;
}

int lcm(int x, int y)
{
  short a, b, c;
  a = x;
  b = y;
  while ((c = a % b) != 0) {
    a = b;
    b = c;
  }
  return (x * y / b);
}

int rp(short * ptr)
{
  short i, j;
  for (i = 1; i <= npt; i++)
    orb[i] = 0;
  for (i = 1; i <= npt; i++) {
    scanf("%hd", ptr + i);
    j = ptr[i];
    if (j <= 0 || j > npt || orb[j])
      return (2);
    orb[j] = 1;
  }
  return (0);
}

int allorbs(void)
{
  short orct, lo, u, v, w, x, y, z;
  for (u = 1; u <= npt; u++)
    orno[u] = 0;
  orct = 0;
  for (u = 1; u <= npt; u++)
    if (orno[u] == 0) {
      orct++;
      orno[u] = orct;
      lo = 1;
      orb[1] = u;
      for (x = 1; x <= lo; x++) {
        z = orb[x];
        for (y = 1; y <= *cp; y++) {
          w = cp[y];
          v = pptr[w][z];
          if (orno[v] == 0) {
            lo++;
            orb[lo] = v;
            orno[v] = orct;
          }
        }
      }
      printf("Orbit number %3d,  length %3d:\n", orct, lo);
      for (x = 1; x <= lo; x++)
        printf("%4d", orb[x]);
      printf("\n");
    }
  return (0);
}

int num_digits(int n)
{
  if (n < 10)
    return (1);
  if (n < 100)
    return (2);
  if (n < 1000)
    return (3);
  if (n < 10000)
    return (4);
  if (n < 100000)
    return (5);
  return (6);
}
