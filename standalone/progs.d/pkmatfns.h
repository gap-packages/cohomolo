#ifndef COHOMOLO_PKMATFNS_H
#define COHOMOLO_PKMATFNS_H

#include "defs.h"

extern char  prime, **mat[], cvec[], pinv[], aut;
extern short dim, svec[], maxnull;
extern FILE *ip, *op;

void trans(char ** a, char ** b);
void copy(char ** a, char ** b);
void ncopy(int n, char ** a, char ** b);
void sum(int n, char ** a, char ** b);
void im(char * v, char * w, char ** a);
void prod(char ** a, char ** b, char ** c);
void readmat(char ** a);
void rvecsum(int n, char * v);
void printmat(char ** a);
int  null(char ** a, char ** b);
int  null1(char ** a, char ** b);
int  spgen(char ** a, int n);
void opnmat(char ** a, int n, int tdim, int fop);

#endif
