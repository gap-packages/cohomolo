#include "defs.h"

#define PSP 2000000
#define SVSP 600000
#define NPT 32767
#define MP 1000
#define MB 80

char heqg, words, fullg, check, inf1[80], inf2[80], inf3[80], outf0[80],
    outf[80];
/* inf1 = gpname.outperm = inf2
   inf3=gpname.inperm (in default case for inf1 only)
   inf1=inf2 if -e set
   outf defined by user of program interactively
   outf0 used to remember groupname.
*/
short perm[PSP], sv[SVSP], cp[10 * NPT], actgen[MP], orb[NPT + 1], base[MB],
    lorbc[MB], lorbg[MB], lorbh[MB], pno[MP / 2], *pptr[MP], *svptr1[MB],
    *svptr2[MB], *svptr3[MB], scpf[MB], scps[MB], sadpt[MB], mp = MP,
                                                             mb = MB - 1;
int psp = PSP, svsp = SVSP;

int main(int argc, char * argv[])
{
  short arg;
  char  c, d, err;
  check = err = words = heqg = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while ((c = argv[arg][0]) == '-') {
    d = argv[arg][1];
    arg++;
    if (argc <= arg) {
      err = 1;
      goto error;
    }
    if (d == 'e')
      heqg = 1;
    else if (d == 'w')
      words = 1;
    else if (d == 't')
      check = 1;
    else {
      err = 1;
      goto error;
    }
  }
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  strcpy(outf0, inf1);
  strcpy(inf2, inf1);
  strcpy(inf3, inf1);
  strcat(inf3, "inperm");
  if (heqg == 0) {
    arg++;
    if (argc <= arg)
      heqg = 1;
    else
      strcat(inf2, argv[arg]);
  }
  arg++;
  if (argc <= arg) {
    fullg = 1;
    strcat(inf1, "outperm");
  }
  else {
    strcat(inf1, argv[arg]);
    fullg = 0;
  }
  if (optprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:   optrun [-e] [-w] [-t] gpname [inf2] [inf1].\n");
    exit(1);
  }
  exit(0);
}
