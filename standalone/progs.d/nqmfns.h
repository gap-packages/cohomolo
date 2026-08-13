#ifndef COHOMOLO_NQMFNS_H
#define COHOMOLO_NQMFNS_H

#include "defs.h"

extern char  inf1[], outf[], outfm[], gap;
extern short intexp, mexp, mng, wksp, prime, expo, nng, class, *rpf, *rpb,
    *eexpnt, *enexpnt, **pcb, mnng, mord, rel[], expnt[], nexpnt[], cord[],
    wt[], d1[], d2[], *pcptr[], **powptr[], **comptr[], *sspc[], *sspf[],
    sgen[], sex[], spgen[], spex[], spugen[], *tlintg[];
extern int ptrsp, rsp;

int  ingp(void);
void outgp(void);
void zero(short * p1, short * p2);
void setnr(short * p);
int  collect(short * spc, short * spf, int sgn);
int  intgen(int i, int j);
int  subrel(int i, int j);
int  assoc(int g1, int g2, int g3);
int  prnrel(int corrtl);

#endif
