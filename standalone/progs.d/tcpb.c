#include "defs.h"

extern char ginrel[];
extern int  space;
extern int  mpt, rel[], cosno[], gno[], inv[], gch[], *imcos[];
int ng, ngi, nsg, nr, endsg, endr, maxcos, str, len, ad, *fpt, *bpt, ccos,
    maxd, totd, lastd, cind, nfree, stcr, endcr, fcos, bcos, lcl;
char   fullsc, clsd, lkah;
FILE * op;
/* All comments as for tcp.c */

void inperr(char * s, int n)
{
  fprintf(stderr, "Input error in %s no %d\n", s, n);
}

int digit(int j)
{
  if (j >= '0' && j <= '9')
    return (1);
  else
    return (0);
}

int letter(int j)
{
  if ((j >= 'a' && j <= 'z') || (j >= 'A' && j <= 'Z'))
    return (1);
  else
    return (0);
}

void seeknln(void)
{
  while (getchar() != '\n')
    ;
}

int readrel(int s, int no)
{
  int  stbr, endbr, exp, l, m, n;
  char ch;
  char gotg, br, clbr, emptybr, name[10];
  if (s)
    strcpy(name, "relation");
  else
    strcpy(name, "subgen");
  gotg = 0;
  br = 0;
  clbr = 0;
  ad = str;
  len = 0;
  ch = getchar();
  while (ch != '\n') {
    if (ch == ' ')
      ch = getchar();
    else if (ch == '(') {
      if (br || clbr) {
        inperr(name, no);
        return (-1);
      }
      emptybr = 1;
      br = 1;
      stbr = ad + 1;
      gotg = 0;
      ch = getchar();
    }
    else if (ch == ')') {
      if (br == 0 || emptybr) {
        inperr(name, no);
        return (-1);
      }
      br = 0;
      gotg = 0;
      clbr = 1;
      endbr = ad;
      ch = getchar();
    }
    else if (letter(ch)) {
      if (clbr || gch[ch] == 0) {
        inperr(name, no);
        return (-1);
      }
      emptybr = 0;
      gotg = 1;
      len++;
      ad++;
      rel[ad] = gch[ch];
      ch = getchar();
    }
    else if (ch == '-' || digit(ch)) {
      if (gotg == 0 && clbr == 0) {
        inperr(name, no);
        return (-1);
      }
      gotg = 0;
      if (ch == '-') {
        ch = getchar();
        if (digit(ch) == 0)
          exp = -1;
        else {
          exp = 0;
          while (digit(ch)) {
            exp *= 10;
            exp -= (ch - '0');
            ch = getchar();
          }
        }
      }
      else {
        exp = 0;
        while (digit(ch)) {
          exp *= 10;
          exp += (ch - '0');
          ch = getchar();
        }
      }
      if (exp == 0) {
        inperr(name, no);
        return (-1);
      }
      if (clbr) {
        if (exp < 0)
          for (m = stbr, n = endbr; m <= n; m++, n--) {
            if (m == n)
              rel[m] = -rel[m];
            else {
              l = -rel[m];
              rel[m] = -rel[n];
              rel[n] = l;
            }
          }
        exp = abs(exp);
        exp--;
        clbr = 0;
        for (n = 1; n <= exp; n++)
          for (m = stbr; m <= endbr; m++) {
            len++;
            ad++;
            rel[ad] = rel[m];
          }
      }
      else {
        n = rel[ad];
        if (exp < 0) {
          n = -n;
          rel[ad] = n;
          exp = -exp;
        }
        exp--;
        for (m = 1; m <= exp; m++) {
          len++;
          ad++;
          rel[ad] = n;
        }
      }
    }
    else {
      inperr(name, no);
      return (-1);
    }
  }
  if (br || clbr) {
    inperr(name, no);
    return (-1);
  }
  rel[str] = len;
  str = ad + 1;
  return (0);
}

int tcprog(void)
{
  int  i, j, k, ch, *p;
  char fname[80];
  printf("Todd-Coxeter Coset Enumeration Algorithm. (HLT + Lookahead).\n\n");
  printf("    Total space for relations and coset table = %d.\n", space);
  printf("    Generators should be upper or lower case letters.\n");
  printf("    Relations and subgroup generators are input as strings in the");
  printf(" generators.\n");
  printf("    Generators may be followed by a positive or negative exponent");
  printf(" ('-' = '-1').\n");
  printf(
      "    Brackets (but no nested brackets) are allowed, and a bracketed");
  printf(" expression.\n");
  printf("      must be followed by an exponent.\n");
  printf("    Example:  ab-3c6(ab)-2(XYZ)-ba-\n\n\n");
  printf("Input numbers of generators, subgroup generators and relations.\n");
  scanf("%d%d%d", &ng, &nsg, &nr);
  for (i = 0; i <= 127; i++)
    gch[i] = 0;
  printf("Input names of generators (upper or lower case letters).\n");
  for (i = 1; i <= ng; i++) {
    ch = ' ';
    while (ch == ' ' || ch == '\n')
      ch = getchar();
    if (letter(ch) == 0) {
      fprintf(stderr, "Generators must be letters.\n");
      return (-1);
    }
    if (gch[ch] != 0) {
      fprintf(stderr, "Repeated generator.\n");
      return (-1);
    }
    gch[ch] = i;
    gno[i] = 0;
    ginrel[i] = 0;
  }
  seeknln();
  if (nsg > 0)
    printf("Input subgroup generators (separated by new lines).\n");
  str = 0;
  ad = -1;
  for (i = 1; i <= nsg; i++)
    if (readrel(0, i) == -1)
      return (-1);
  endsg = ad;
  if (nr > 0)
    printf("Input relators (separated by new lines).\n");
  for (i = 1; i <= nr; i++) {
    if (readrel(1, i) == -1)
      return (-1);
    if (len == 2 && rel[ad] == rel[ad - 1]) {
      gno[abs(rel[ad])] = 1;
      str -= 3;
      ad -= 3;
    }
    else {
      k = ad - len + 1;
      for (j = k; j <= ad; j++)
        ginrel[abs(rel[j])] = 1;
    }
  }
  for (i = 1; i <= ng; i++)
    if (gno[i] == 1 && ginrel[i] == 0) {
      rel[ad + 1] = 2;
      rel[ad + 2] = i;
      rel[ad + 3] = i;
      ad += 3;
      gno[i] = 0;
    }
  endr = ad;

  ngi = -1;
  for (i = 1; i <= ng; i++)
    if (gno[i] == 1) {
      ngi++;
      gno[i] = ngi;
      gno[i + 52] = ngi;
      inv[ngi] = ngi;
    }
    else {
      ngi += 2;
      gno[i] = ngi - 1;
      gno[i + 52] = ngi;
      inv[ngi] = ngi - 1;
      inv[ngi - 1] = ngi;
    }
  i = -1;
  while (++i <= endr) {
    j = rel[i];
    k = i + j;
    while (++i <= k) {
      p = rel + i;
      *p = (*p < 0) ? gno[52 - *p] : gno[*p];
    }
    i = k;
  }
  maxcos = (space - endr - 2) / (ngi + 3);
  fpt = rel + 1 + endr;
  bpt = fpt + maxcos;
  for (i = 0; i <= ngi; i++)
    imcos[i] = fpt + maxcos * (2 + i);
  printf("Maxcos=%d.\n\n", maxcos);

  ccos = 1;
  lastd = 1;
  maxd = 1;
  totd = 1;
  cind = 1;
  nfree = 2;
  for (i = 0; i <= ngi; i++)
    imcos[i][1] = 0;
  for (i = 0; i < maxcos; i++)
    fpt[i] = i + 1;
  fpt[1] = 0;
  fpt[maxcos] = 0;
  lkah = 0;
  bpt[1] = 0;
  endcr = -1;
  while (endcr != endsg) {
    scanrel();
    if (fullsc == 0) {
      fprintf(stderr, "Space overflow in subgroup generation phase.\n");
      return (-1);
    }
  }
  while (ccos != 0) {
    clsd = 1;
    endcr = endsg;
    while (endcr != endr) {
      scanrel();
      if (fullsc == 0) {
        clsd = 0;
        if (lkah == 0) {
          lcl = bpt[ccos];
          printf("Entering lookahead.\n");
          lkah = 1;
        }
      }
    }
    if (lkah) {
      i = fpt[ccos];
      if (clsd) {
        j = bpt[ccos];
        if (j != lcl) {
          fpt[j] = i;
          if (i == 0)
            lastd = j;
          else
            bpt[i] = j;
          j = fpt[lcl];
          fpt[lcl] = ccos;
          bpt[ccos] = lcl;
          fpt[ccos] = j;
          bpt[j] = ccos;
        }
        lcl = ccos;
      }
      ccos = i;
      if (ccos == 0) {
        if (cind == maxcos) {
          fprintf(stderr, "Not enough space.\n");
          return (-1);
        }
        printf("Exiting lookahead. No. of cosets=%d\n", cind);
        ccos = fpt[lcl];
        lkah = 0;
      }
    }
    else
      ccos = fpt[ccos];
  }

  for (i = 1; i <= ng; i++)
    if (ginrel[i] == 0) {
      k = gno[i];
      j = lastd;
      while (j != 0) {
        if (imcos[k][j] == 0) {
          fprintf(
              stderr,
              "Coset table incomplete at end of scan. Index is infinite.\n");
          return (-1);
        }
        j = bpt[j];
      }
    }
  printf("Algorithm complete. maxdef,totdef=%d,%d.\n\n", maxd, totd);
  if (nsg == 0)
    printf("The order of the group is %d.\n\n", cind);
  else
    printf("The index of the subgroup is %d.\n\n", cind);
  if (cind <= mpt) {
    printf("Do you wish to store the permutations? (y/n)    ");
    ch = getchar();
    seeknln();
    if (ch == 'y') {
      printf("Input filename    ");
      scanf("%s", fname);
      op = fopen(fname, "w");
      fprintf(op, "%4d%4d%4d%4d\n", cind, ng, 0, 0);
      bpt[1] = 1;
      cosno[1] = 1;
      j = 1;
      for (i = 2; i <= cind; i++) {
        j = fpt[j];
        cosno[i] = j;
        bpt[j] = i;
      }
      for (i = 1; i <= ng; i++) {
        k = gno[i];
        for (j = 1; j <= cind; j++)
          fprintf(op, " %3d", bpt[imcos[k][cosno[j]]]);
        fprintf(op, "\n");
      }
    }
  }
  return (0);
}

int scanrel(void)
{
  int  i, j, k, l, m;
  char comp;
  fullsc = 1;
  stcr = endcr + 2;
  endcr += (1 + rel[stcr - 1]);
  fcos = ccos;
  bcos = ccos;
  comp = 1;
  for (i = stcr; i <= endcr; i++) {
    k = imcos[rel[i]][fcos];
    if (k == 0) {
      comp = 0;
      break;
    }
    fcos = k;
  }
  if (comp) {
    if (fcos != bcos)
      coinc(fcos, bcos);
    return (0);
  }
  stcr = i;
  for (i = endcr; i >= stcr; i--) {
    l = rel[i];
    m = inv[l];
    k = imcos[m][bcos];
    if (k == 0) {
      if (i == stcr) {
        k = imcos[l][fcos];
        if (k == 0) {
          imcos[l][fcos] = bcos;
          imcos[m][bcos] = fcos;
        }
        else if (k != bcos)
          coinc(k, bcos);
        return (0);
      }
      if (lkah || nfree == 0) {
        fullsc = 0;
        return (0);
      }
      for (j = 0; j <= ngi; j++)
        imcos[j][nfree] = 0;
      totd++;
      cind++;
      if (maxd < cind)
        maxd++;
      imcos[m][bcos] = nfree;
      imcos[l][nfree] = bcos;
      bcos = nfree;
      bpt[nfree] = lastd;
      fpt[lastd] = nfree;
      lastd = nfree;
      nfree = fpt[nfree];
      fpt[lastd] = 0;
    }
    else
      bcos = k;
  }
  if (fcos != bcos)
    coinc(fcos, bcos);
  return (0);
}

int coinc(int c1, int c2)
{
  int lc, hc, qh, qt, i, j, x, fhc, bhc, lim, him;
  if (c1 < c2) {
    lc = c1;
    hc = c2;
  }
  else {
    lc = c2;
    hc = c1;
  }
  qh = 0;
  qt = 0;
  fhc = fpt[hc];
  bhc = bpt[hc];
  fpt[bhc] = fhc;
  if (fhc == 0)
    lastd = bhc;
  else
    bpt[fhc] = bhc;
  if (ccos == hc) {
    ccos = bhc;
    endcr = endr;
    clsd = 0;
  }
  if (lkah && lcl == hc)
    lcl = bhc;
  while (1) {
    fpt[hc] = nfree;
    nfree = hc;
    cind--;
    for (i = 0; i <= ngi; i++) {
      him = imcos[i][hc];
      if (him != 0) {
        j = inv[i];
        lim = imcos[i][lc];
        if (him == hc)
          him = lc;
        else
          imcos[j][him] = 0;
        if (lim == 0)
          imcos[i][lc] = him;
        else {
          if (lim == hc) {
            imcos[i][lc] = lc;
            lim = lc;
          }
          while (fpt[him] < 0)
            him = -fpt[him];
          while (fpt[lim] < 0)
            lim = -fpt[lim];
          if (him != lim) {
            if (lim > him) {
              x = lim;
              lim = him;
              him = x;
            }
            fhc = fpt[him];
            bhc = bpt[him];
            fpt[bhc] = fhc;
            if (fhc == 0)
              lastd = bhc;
            else
              bpt[fhc] = bhc;
            if (ccos == him) {
              ccos = bhc;
              endcr = endr;
              clsd = 0;
            }
            if (lkah && lcl == him)
              lcl = bhc;
            fpt[him] = -lim;
            if (qh == 0)
              qh = him;
            else
              bpt[qt] = him;
            qt = him;
            bpt[qt] = 0;
          }
        }
        x = imcos[i][lc];
        if (imcos[j][x] == 0)
          imcos[j][x] = lc;
      }
    }
    if (qh == 0)
      break;
    hc = qh;
    qh = bpt[qh];
    lc = -fpt[hc];
    bpt[hc] = 0;
  }
  return (0);
}
