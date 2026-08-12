#ifndef COHOMOLO_EXC_H
#define COHOMOLO_EXC_H

int gc(void);
int comp(short * a, short ** dp);
int compb(short * a, short ** dp);
int expand(short * a, short * b);
int action(short * a, short * b);
int concat(short * a, short * b);
int concatl(short * b);
int invwd(short * a, short * b);
int ainvb(short * a, short * b, short * c);
int reduce(short * a);
int concheck(void);
int setpinv(void);

#endif
