#ifndef COHOMOLO_CHB_H
#define COHOMOLO_CHB_H

#include "defs.h"

extern short npt, nf, cp[], orb[], pno[], fp[], tsv1[], tsv2[], tsv3[],
    orep[], *pptr[], mb;

int addperm(void);
int intbase(int      pt,
            int      pos,
            short *  stad,
            short *  nbad,
            short *  b,
            short *  lorb,
            short ** svptr);

#endif
