#ifndef COHOMOLO_MCP_H
#define COHOMOLO_MCP_H

void seeknln(void);
void addsv(int pt, short * sv);
void invert(short * ptr1, short * ptr2);
int  image(int pt);
void readvec(short * ptr, int e);
void readbaselo(int nb, short * base, short * lorb);
void readpsv(int e, int nb, int nperms, short ** svptr);
void setpinv(void);
int  mcprog(void);
int  conprog(int con);

#endif
