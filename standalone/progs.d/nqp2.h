#ifndef COHOMOLO_NQP2_H
#define COHOMOLO_NQP2_H

void enter(short * g, int pow);
void entvec(short * h, short * g, int pow);
void expand(short * p1, short * p2, int len);
void compress(short * p1, short * p2, int len);
void nchg(void);
int  spact(void);
int  intact(void);

#endif
