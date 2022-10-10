#include "defs.h"

#define SPACE 15000000
#define MPT 4000
/* All comments as for tcd.c  */

char rs, ginrel[53];
int  space = SPACE;
int  mpt = MPT, rel[SPACE], cosno[MPT + 1], gno[105], inv[104], gch[128],
    *imcos[104];

int main(void)
{
  rs = 0;
  if (tcprog() == -1)
    exit(1);
  exit(0);
}
