#ifndef COHOMOLO_MOREPERMFNS_H
#define COHOMOLO_MOREPERMFNS_H

#include "defs.h"

extern short npt, pno[], expcp[], cp[], genorb[], *pptr[], *expptr[];

int expandp(int nb, short * p, short * base, short ** svptr);
int resetsv(int nb, short * base, short * lorb, short * gno, short ** svptr);
int allorbs(short * lorb, short * orno);
int backimage(int pt);
int exprep(int pt, int no, short * sv);
int expimage(int pt);

#endif
