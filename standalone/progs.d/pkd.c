#include "defs.h"

#define SPACE 6000000
#define MDIM 2000
#define MV 20000
#define MM 10
#define MPR 127
#define MSL 50

/* MPR and MM cannot sensibly be increased. Matrices are stored as
   characters, so entries must lie in range -128<=x<=127. MM=10, means
   that groups can have at most 9 generators. This is only to allow the
   gen no. to be read as a single digit.
*/
char inf[80], outf[80], temp[80], temp2[80], mspace[SPACE], *vec[MV],
    **mat[MM],
    /* defaults inf = gpname.inmat
       outf always takes values infq infs or inff.
    */
    cvec[MDIM + 1], pinv[MPR], mstr[MSL + 1], full, intop, opt, aut;
short svec[MDIM + 1], mdim = MDIM, mv = MV, mm = MM, mpr = MPR, msl = MSL;
int   space = SPACE;

int main(int argc, char * argv[])
{
  char arg, err, c;
  opt = 0;
  intop = 0;
  full = 0;
  aut = 0;
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    if ((c = argv[arg][1]) == 'f')
      full = 1;
    else if (c == 'i')
      intop = 1;
    else if (c == 'o')
      opt = 1;
    else if (c == 'a')
      aut = 1;
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
  strcpy(temp, inf);
  strcat(temp, "temp");
  strcpy(temp2, temp);
  strcat(temp2, "2");
  arg++;
  if (argc <= arg)
    strcat(inf, "inmat");
  else
    strcat(inf, argv[arg]);
  strcpy(outf, inf);
  if (pkprog() == -1)
    exit(1);
error:
  if (err) {
    printf("Usage: pkrun [-f] [-i] [-o] gpname [inf]\n");
    exit(1);
  }
  exit(0);
}
