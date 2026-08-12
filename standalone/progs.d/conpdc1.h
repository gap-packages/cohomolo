#ifndef COHOMOLO_CONPDC1_H
#define COHOMOLO_CONPDC1_H

int  image(int pt);
void addsvb(int pt, int ** sv);
void addsvf(int pt, int ** sv);
void invert(int * a, int * b);
void rdperm(int * a, int * b);
void rdsv(int ** sv);
int  cnprg1(void);
void advance(void);

#endif
