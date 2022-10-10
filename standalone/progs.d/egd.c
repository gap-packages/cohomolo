#include "defs.h"

#define NPT 32767
#define MXP 501
#define MNB 80
/* MXP-1 = max no of perms. */

char inf[80], inf2[80], outf[80], full, stabcall;
/* defaults inf gpname.outperm
           outf gpname.sg
*/
short perm[(MXP - 1) * (NPT + 1)], sv2[(MNB - 1) * NPT], cp[1], orb[NPT + 1],
    base[MNB], lorb[MNB], fixpt[MNB + 1], fixb[MXP], pno[MXP], *pptr[MXP],
    *svptr[MNB], mnpt = NPT, mp = MXP - 1, mb = MNB - 1;

int main(int argc, char * argv[])
{
  short arg;
  char  err, c;
  stabcall = err = 0;
  arg = 1;
  full = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    if ((c = argv[arg][1]) == 'f')
      full = 1;
    else if (c == 's')
      stabcall = 1;
    else {
      err = 1;
      goto error;
    }
    arg++;
    if (argc <= arg) {
      err = 1;
      goto error;
    }
  }
  strcpy(inf, argv[arg]);
  strcat(inf, ".");
  strcpy(outf, inf);
  strcpy(inf2, inf);
  arg++;
  if (argc <= arg)
    strcat(inf, "outperm");
  else
    strcat(inf, argv[arg]);
  if (full) {
    arg++;
    if (argc <= arg)
      strcat(inf2, "inperm");
    else
      strcat(inf2, argv[arg]);
  }
  arg++;
  if (argc <= arg)
    strcat(outf, "sg");
  else
    strcat(outf, argv[arg]);
  if (egprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(
        stderr,
        "Usage:    egrun [-s] [-f] gpname [inf] (if -f [inf2]) [outf].\n");
    exit(1);
  }
  exit(0);
}
