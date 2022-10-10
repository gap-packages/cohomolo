#include "defs.h"

#define NPT 2000
int   npt, perm[NPT + 1], had[NPT + 1];
FILE *fopen(), *op;

int main(int argc, char * argv[])
{
  int  arg, nb, nperms, i, j, k, l, n;
  char cyc, outf[80], err;
  int  c;
  /* Default outf=gpname.inperm  */
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  strcpy(outf, argv[arg]);
  strcat(outf, ".");
  arg++;
  if (argc <= arg)
    strcat(outf, "inperm");
  else
    strcat(outf, argv[arg]);
  op = fopen(outf, "w");
  printf("Input npt, nperms, init. nb.\n");
  scanf("%d%d%d", &npt, &nperms, &nb);
  if (npt > NPT) {
    fprintf(stderr, "npt too big. Increase NPT.\n");
    exit(1);
  }
  fprintf(op, "%4d%4d%4d%4d\n", npt, nperms, nb, 0);
  if (nb != 0) {
    printf("Input initial base points.\n");
    for (i = 1; i <= npt; i++)
      had[i] = 0;
    for (i = 1; i <= nb; i++) {
      scanf("%d", &n);
      if (n <= 0 || n > npt || had[n]) {
        fprintf(stderr, "Invalid or repeated base point %d\n", n);
        exit(1);
      }
      perm[i] = n;
      had[n] = 1;
    }
    if (npt >= 1000) {
      for (i = 1; i <= nb; i++)
        fprintf(op, "%5d", perm[i]);
      fprintf(op, "\n");
    }
    else {
      for (i = 1; i <= nb; i++)
        fprintf(op, "%4d", perm[i]);
      fprintf(op, "\n");
    }
  }
  printf("Now input perms in cyclic not'n. End each perm with a '.'.\n");
  for (i = 1; i <= nperms; i++) {
    for (j = 1; j <= npt; j++) {
      perm[j] = j;
      had[j] = 0;
    }
    while ((c = getchar()) != '.')
      if (c == '(') {
        cyc = 1;
        j = 0;
        while (cyc) {
          scanf("%d", &l);
          if (l <= 0 || l > npt || had[l]) {
            fprintf(stderr, "Invalid or repeated point %d\n", l);
            exit(1);
          }
          had[l] = 1;
          if (j == 0) {
            j = l;
            k = l;
          }
          else {
            perm[k] = l;
            k = l;
          }
          while ((c = getchar()) == ' ')
            ;
          if (c == ')') {
            cyc = 0;
            perm[k] = j;
          }
        }
      }
    if (npt >= 1000) {
      for (j = 1; j <= npt; j++)
        fprintf(op, "%5d", perm[j]);
      fprintf(op, "\n");
    }
    else {
      for (j = 1; j <= npt; j++)
        fprintf(op, "%4d", perm[j]);
      fprintf(op, "\n");
    }
  }
error:
  if (err) {
    fprintf(stderr, "Usage    makegp gpname [outf].\n");
    exit(1);
  }
  exit(0);
}
