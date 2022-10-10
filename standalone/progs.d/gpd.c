#include "defs.h"

#define NPT 32767
#define PSP 2000000
#define SVSP 200000
#define MP 1000
#define MB 80

/* NPT = max no of points,  MP = max no of perms,  PSP = space for perms,
   SVSP = space for Schreier vectors,  MB-1 = max no of base pts.
*/
char wrd, nt, isbase, inf[80], outf1[80], outf2[80], fixed[NPT + 1];
/* defaults: inf gpname.inperm
             outf1 gpname.outperm
             outf2 gpname.words
*/
short perm[PSP], sv[SVSP], cp[5 * NPT], actgen[MP], orb[NPT + 1], base[MB],
    lorb[MB], order[MB], pno[MP / 2], *pptr[MP], *svptr[MB],
    mp = MP, mb = MB - 1, mnpt = NPT;
int psp = PSP, svsp = SVSP;

int main(int argc, char * argv[])
{
  short arg;
  char  d, err;
  nt = wrd = isbase = 0;
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    d = argv[arg][1];
    arg++;
    if (d == 'w')
      wrd = 1;
    else if (d == 'n')
      nt = 1;
    else if (d == 'b')
      isbase = 1;
    else {
      err = 1;
      goto error;
    }
    if (argc <= arg) {
      err = 1;
      goto error;
    }
  }
  strcpy(inf, argv[arg]);
  strcat(inf, ".");
  strcpy(outf1, inf);
  if (wrd) {
    strcpy(outf2, inf);
    strcat(outf2, "words");
  }
  arg++;
  if (argc <= arg)
    strcat(inf, "inperm");
  else
    strcat(inf, argv[arg]);
  arg++;
  if (argc <= arg)
    strcat(outf1, "outperm");
  else
    strcat(outf1, argv[arg]);
  if (gpprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:   gprun [-n] [-w] [-b] gpname [inf] [outf1].\n");
    exit(1);
  }
  else
    exit(0);
}
