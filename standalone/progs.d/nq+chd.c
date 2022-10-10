#include "defs.h"

#define MSP 60000
#define MM 100
#define MV 2000
#define MDIM 100
#define MPR 256
#define MWL 200
int   msp = MSP;
short mv = MV, mm = MM, mdim = MDIM, mpr = MPR, prime, dim, *spv, **spm,
      ngens, mspace[MSP], *vec[MV], **mat[MM], cp[MWL], pinv[MPR], opmats;
char   inf1[80], inf2[80];
FILE * ip;

int main(int argc, char * argv[])
{
  short arg;
  char  err;
  err = 0;
  arg = 1;
  opmats = 0;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  if (argv[arg][0] == '-') {
    if (argv[arg][1] == 'o')
      opmats = 1;
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
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  strcpy(inf2, inf1);
  if (calcmats() == -1)
    exit(1);
  if (matact(0) == -1)
    exit(1);

  arg++;
  while (argc > arg) {
    strcpy(inf1, inf2);
    strcat(inf1, argv[arg]);
    if (matact(1) == -1)
      exit(1);
    arg++;
  }
error:
  if (err) {
    fprintf(stderr, "Usage:  nq+chrun gpname [inf] [inf] ...\n");
    exit(1);
  }
  exit(0);
}
