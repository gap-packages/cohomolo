#ifndef COHOMOLO_PERMFNS_H
#define COHOMOLO_PERMFNS_H

short orbitsv(short pt, short * sv, short lo);
short addsv(short pt, short * sv);
short image(short pt);
short invert(short * ptr1, short * ptr2);
short readperm(short * ptr);
short printvec(short * ptr, short e);
short readvec(short * ptr, short e);
short readbaselo(short nb, short * base, short * lorb);
short printbaselo(short nb, short * base, short * lorb);
short printpsv(short nb, short * gno, short ** svptr);
short readpsv(short e, short nb, short nperms, short ** svptr);

#endif
