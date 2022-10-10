#include "defs.h"

#define PSP 500000
#define SVSP 100000
#define NPT 32767
#define MB 80
#define TRSP 2000000
/* TRSP is the space of the array used to store the coset reps, which are
   stored using a tree structure.
*/

char inf0[80], inf1[80], inf2[80], inf3[80], outf1[80], outf2[80], outf3[80],
    temp1[90], temp2[90],
    /*  inf1 is gpname.?1 (no default)
        inf2 is gpname.?2 (no default, the subgroup of inf1 for which
       transversal is to be computed). Defaults: inf3  gpname.pg outf1 inf3.nr
                 outf2  gpname.dcr
                 outf3  gpname.cr.
        inf0 is used to remember the gpname alone.
    */
    ofstr[2], ngn[8], expg, exph, cr, dcr, hg, triv, cong;
int psp = PSP, trsp = TRSP, svsp = SVSP;
int mpt = NPT, mb = MB, nfuse, tree[TRSP], perm[PSP], gorb[NPT + 1],
    lorbg[MB], lorbh[MB], base[MB], fpt[MB], bpt[MB], coh_index[MB + 1],
    *cp[10 * NPT], *trad[MB], *sv[SVSP], *tailad[MB], **cpad[MB],
    **svgptr[MB], **svhptr[MB];

int main(int argc, char * argv[])
{
  int  arg;
  char c, err;
  int  d;
  err = 0;
  cong = 0;
  arg = 1;
  hg = 0;
  cr = 0;
  expg = 0;
  exph = 0;
  dcr = 0;
  nfuse = -1;
  if (argc <= arg) {
    err = 1;
    goto error;
  }
  while (argv[arg][0] == '-') {
    c = argv[arg][1];
    if (c == 'c')
      cr = 1;
    else if (c == 'g')
      cong = 1;
    else if (c == 'e') {
      d = argv[arg][2];
      if (d == 'g')
        expg = 1;
      else if (d == 'h')
        exph = 1;
      else {
        err = 1;
        goto error;
      }
    }
    else if (c == 'd') {
      dcr = 1;
      hg = 1;
      d = argv[arg][2];
      if (d != '\0') {
        d -= '0';
        if (d < 0 || d > 9) {
          err = 1;
          goto error;
        }
        ofstr[1] = '\0';
        nfuse = d;
      }
    }
    else if (c == 'h')
      hg = 2;
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
  if (argc <= arg + 2) {
    err = 1;
    goto error;
  }
  strcpy(inf1, argv[arg]);
  strcat(inf1, ".");
  strcpy(inf2, inf1);
  strcpy(inf3, inf1);
  strcpy(outf1, inf1);
  strcpy(outf3, inf1);
  strcpy(outf2, inf1);
  strcpy(inf0, inf1);
  strcpy(temp1, inf1);
  strcpy(temp2, inf1);
  strcat(temp1, "conpt");
  strcat(temp2, "crt");
  strcat(inf1, argv[arg + 1]);
  arg += 2;
  if (strcmp(argv[arg], "triv") == 0)
    triv = 1;
  else {
    triv = 0;
    strcat(inf2, argv[arg]);
  }
  if (hg) {
    arg++;
    if (argc <= arg)
      strcat(inf3, "pg");
    else
      strcat(inf3, argv[arg]);
  }
  if (dcr) {
    arg++;
    if (argc <= arg)
      strcat(outf2, "dcr");
    else
      strcat(outf2, argv[arg]);
  }
  if (cr) {
    arg++;
    if (argc <= arg)
      strcat(outf3, "cr");
    else
      strcat(outf3, argv[arg]);
  }
  if (dcr && nfuse >= 0) {
    arg++;
    if (argc <= arg)
      strcpy(ngn, "ng");
    else
      strcpy(ngn, argv[arg]);
  }
  if (hg && dcr == 0) {
    arg++;
    if (argc <= arg) {
      strcpy(outf1, inf3);
      strcat(outf1, ".nr");
    }
    else
      strcat(outf1, argv[arg]);
  }
  if (cnprg1() == -1)
    exit(1);
  if (hg || cong)
    if (cnprg2() == -1)
      exit(1);
error:
  if (err) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "conrun [-g] [-c] [-d(n)] [-h] [-eg] [-eh]\n");
    fprintf(stderr,
            "gpname inf1 inf2 (if -h [inf3]) [outf2(dcr)] [outf3(cr)].\n");
    exit(1);
  }
  exit(0);
}
