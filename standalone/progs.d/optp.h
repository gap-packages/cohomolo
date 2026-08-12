#ifndef COHOMOLO_OPTP_H
#define COHOMOLO_OPTP_H

int  optprog(void);
int  gp(int spno, int sbno, short * lorb, short ** svptr);
int  test(short ** svp, int expo);
int  nc(int th, int sno, int eno);
int  comm(void);
int  outp(int c);
int  intsect(short ** sv1, short ** sv2, short ** sv3, int stint);
int  im(int pt);
void addsvb(int pt, short * sv);
void addsvf(int pt, short * sv);
int  core(void);
void rhc(void);
void rgh(int c);
int  join(void);
int  spacecheck(void);

#endif
