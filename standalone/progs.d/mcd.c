#include "defs.h"

#define MSP 2500000
#define MM 400
#define MV 25000
#define NPT 4000
#define PSP 1000000
#define SVSP 75000
#define MP 500
#define MPR 500
#define MDIM 300
#define MB 80
char slg, con, check, inf1[80], inf2[80], inf3[80], inf4[80], inf5[80],
    inf6[80], outf1[80];
/* inf1 always gpname.inmat
   inf2 always gpname.words
   inf3 gpname.sg (but gpname.outperm if -o set)
   inf4 given as arguments (see Info.3)
Warning: Sequence of arguments must be as described in Info.3, or
         chaos will result.
   inf5 used to remember gpname
   inf6 =gpname.sg.rel if -t set. later =gpname.cbmats.
   outf1 = inf4"mat"
*/
int   msp = MSP, psp = PSP, svsp = SVSP;
short mm = MM, mv = MV, mp = MP, mpt = NPT, mpr = MPR, mdim = MDIM, mb = MB,
      mspace[MSP], *vec[MV], **mat[MM], pinv[MPR], perm[PSP], sv[SVSP],
      cp[5 * NPT], *pptr[MP], *svptr[MB], base[MB], lorb[MB];

int main(int argc, char * argv[])
{
  short arg;
  char  c, err;
  check = 0;
  err = 0;
  arg = 1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  slg = 1;
  c = argv[arg][0];
  while (c == '-') {
    c = argv[arg][1];
    if (c == 'o')
      slg = 0;
    else if (c == 't')
      check = 1;
    else {
      err = 1;
      goto error;
    }
    arg++;
    if (argc <= arg) {
      err = 1;
      goto error;
    }
    c = argv[arg][0];
  }
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  strcpy(inf2, inf1);
  strcpy(inf3, inf1);
  strcpy(inf4, inf1);
  strcpy(inf5, inf1);
  strcpy(inf6, inf1);
  strcat(inf1, "inmat");
  strcat(inf2, "words");
  if (slg)
    strcat(inf3, "sg");
  else
    strcat(inf3, "outperm");
  if (check)
    strcat(inf6, "sg.rel");
  if (mcprog() == -1)
    exit(1);
  while (1)
  /* Generate filenames for conprog
     con is set 1 when argument "pg" is encountered.
     At this stage, all matrices are inverted and transposed, and
     basis changes to the module are made to triangularize the action
     of P, and provide definitions in the PCP of P.d(M)
     if the filename is preceded by a '-', then con is put =2,
     and inverses of matrices are output.
  */
  {
    arg++;
    if (argc <= arg)
      break;
    if (strcmp(argv[arg], "pg") == 0) {
      con = 1;
      strcpy(inf6, inf5);
      strcat(inf6, "cbmats");
    }
    else if (argv[arg][0] == '-')
      con = 2;
    else
      con = 0;
    strcpy(inf4, inf5);
    if (con == 2)
      strcat(inf4, argv[arg] + 1);
    else
      strcat(inf4, argv[arg]);
    strcpy(outf1, inf4);
    strcat(outf1, "mat");
    if (conprog(con) == -1)
      exit(1);
  }
error:
  if (err) {
    fprintf(stderr, "Usage:    matcalc [-o] [-t] gpname\n");
    fprintf(stderr, " [inf4] [inf4] ...\n");
    exit(1);
  }
  else
    exit(0);
}
