#include "defs.h"

#define SPACE 2500000
#define MPT 10000
#define MXC 32383
/* SPACE is total space in coset table.
   MPT is maximum index for which permutations can be output.
   MXC = 2^15-1 is largest positive value of a short integer. Use
   tcrunb if larger values required!
*/

char  rs, ginrel[53];
int   space = SPACE, mxc = MXC;
short mpt = MPT, rel[SPACE], cosno[MPT + 1], gno[105], inv[104], gch[128],
      *imcos[104];

int main(int argc, char * argv[])
{
  char err, arg;
  rs = 0;
  err = 0;
  arg = 1;
  if (argc > arg) {
    if (argv[arg][0] != '-' || argv[arg][1] != 'r') {
      err = 1;
      goto error;
    }
    rs = 1;
    /* -r option sets rs = 1. This allows the user to reduce the amount of
       space in the coset table. This is useful for investigating the effect
       of using lookahead.
    */
  }
  if (tcprog() == -1)
    exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:    tcrun [-r]\n");
    exit(1);
  }
  exit(0);
}
