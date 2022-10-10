#include "defs.h"

#define MP 100
FILE *ip, *op;

void seeknln(void)
{
  while (getc(ip) != '\n')
    ;
}

int main(int argc, char * argv[])
{
  short npt, nb, np, i, j, k, ppt, base, npop, pn, arg;
  char  err, wt, inf[80], outf[80], pno[MP];
  err = 0;
  wt = 0;
  npop = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  if (argv[arg][0] == '-') {
    if (argv[arg][1] == 'w')
      wt = 1;
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
  strcpy(inf, argv[arg]);
  strcat(inf, ".");
  strcpy(outf, inf);
  arg++;
  if (argc <= arg)
    strcat(inf, "spc");
  else
    strcat(inf, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf, "pg");
  else
    strcat(outf, argv[arg]);
  if ((ip = fopen(inf, "r")) == 0) {
    fprintf(stderr, "Cannot open %s.\n", inf);
    goto error;
  }
  op = fopen(outf, "w");
  fscanf(ip, "%hd%hd%hd%hd", &npt, &np, &nb, &k);
  seeknln();
  if (np >= MP) {
    fprintf(stderr, "Increase MP\n");
    goto error;
  }
  if (wt == 0) {
    printf("Input pnos, preceded by their number.\n");
    scanf("%hd", &npop);
  }
  for (i = 1; i <= np; i++)
    pno[i] = 0;
  if (wt) {
    for (i = 1; i <= 2 + np + nb; i++)
      seeknln();
    for (i = 0; i < np; i++) {
      fscanf(ip, "%hd", &j);
      if (j == 1) {
        npop++;
        pno[np - i] = 1;
      }
    }
    fclose(ip);
    ip = fopen(inf, "r");
    seeknln();
  }
  else
    for (i = 1; i <= npop; i++) {
      scanf("%hd", &pn);
      pno[pn] = 1;
    }
  fprintf(op, "%4d%4d%4d%4d\n", npt, npop, nb, 0);
  if (nb > 0) {
    if (npt >= 1000)
      for (i = 1; i <= nb; i++) {
        fscanf(ip, "%hd", &base);
        fprintf(op, "%5d", base);
      }
    else
      for (i = 1; i <= nb; i++) {
        fscanf(ip, "%hd", &base);
        fprintf(op, "%4d", base);
      }
    fprintf(op, "\n");
    seeknln();
    if (k != 0)
      seeknln();
  }
  for (i = 1; i <= np; i++) {
    if (pno[i]) {
      if (npt >= 1000)
        for (j = 1; j <= npt; j++) {
          fscanf(ip, "%hd", &ppt);
          fprintf(op, "%5d", ppt);
        }
      else
        for (j = 1; j <= npt; j++) {
          fscanf(ip, "%hd", &ppt);
          fprintf(op, "%4d", ppt);
        }
      fprintf(op, "\n");
    }
    seeknln();
  }
error:
  if (err) {
    fprintf(stderr, "Usage:    selgen [-w] gpname [inf] [outf] [inf2].\n");
    exit(1);
  }
  else
    exit(0);
}
