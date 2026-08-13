#ifndef COHOMOLO_PCSCFNS_H
#define COHOMOLO_PCSCFNS_H

#include "defs.h"

extern short npt, nb, npt1, expo, prime, pinv[], cp[], power[], wt[], base[],
    ngno[], igno[], *svptr[], *pptr[], d1[];

void firstgen(short * p, short * hg, short * co);
int  express(short * p, short * relc, int nwt);
void setpinv(void);

#endif
