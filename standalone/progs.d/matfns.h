#ifndef COHOMOLO_MATFNS_H
#define COHOMOLO_MATFNS_H

#include "defs.h"

extern short prime, dim, *spv, **spm, **mat[], pinv[];
extern FILE *ip, *op;

void trans(short ** a, short ** b);
void im(short * v, short * w, short ** a);
int  comm(short * v, short * w, short ** a);
void prod(short * cm, short ** a);
int  inv(short ** a, short ** b);
void readmat(short ** a);
void printmat(short ** a);
int  cbdef(int     gb,
           int     ge,
           int     cbno,
           short * d1,
           short * d2,
           short * wt,
           short * acl);

#endif
