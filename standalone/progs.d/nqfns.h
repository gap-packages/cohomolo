#ifndef COHOMOLO_NQFNS_H
#define COHOMOLO_NQFNS_H

#include "defs.h"

extern char  inf[], inf1[], outf1[];
extern short facexp, tails, stage, depth, no, mng, mexp, prime, expo, nng,
    class, dim, onng, *rpf, *rpb, *eexpnt, *enexpnt, **pcb, **opcb, **npcb,
    **npcb2, *nd1, *nd2, **extno, **subno, chpdim, chsdim, rel[], expnt[],
    nexpnt[], prvec[], pinv[], wt[], d1[], d2[], *pcptr[], **powptr[],
    **comptr[], *sspc[], *sspf[], sgen[], sex[], spgen[], spex[], spugen[],
    dpth[];
extern int   rsp, wsp, ptrsp, marg;
extern FILE *ip, *op;

int  ingp(int inp);
int  outgp(void);
void zero(short * p1, short * p2);
void setnr(short * p);
int  collect(short * spc, short * spf, int sgn);
void setpinv(void);
void bgc(void);
int  intgen(int i, int j);
int  subrel(int i, int j);
int  assoc(int g1, int g2, int g3);
int  prnrel(void);

#endif
