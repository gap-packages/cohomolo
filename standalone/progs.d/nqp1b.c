#include "defs.h"

extern char inf0[], inf1[], inf2[], inf3[], inf4[], outf0[], outf1[], outf2[],
    outfd[], inf[], act, ch1, crel, cfm, gap;
extern int mv, mm, facexp, tails, stage, depth, no, mng, mcl, prime, exp, nng,
    class, *rpf, *rpb, **pcb, dim, onng, *spv, **spm, rel[], wt[], d1[], d2[],
    *pcptr[], sd1[], sd2[], swt[], dpth[], mspace[], *vec[], **mat[], cp[];
extern int rsp, msp, ptrsp, wsp;
int **     intg, **cintg, cbno, ngens, maxm, maxv, matcl, **extno, **subno,
    chsdim = 0, chpdim, exp1;
int   orsp, optrsp, rspk;
long  inf3offset, inf4offset;
char  norm;
FILE *ip, *ipm, *op;

/* General comments on programs nqrun (and nqmrun).
   See file Info.5 for format of ip/op files, and meaning of variables
   exp,prime,facexp,no,class,depth.
   Definiton of i is [d1[i],d2[i]] or d1[i]^prime if equal.
   wt[i] and dpth[i] are weights and depths.
   sd1,sd2 and swt are used to recall defs and weights of the calculation
   in a Sylow intersection.
   rel is array used to hold all pcp definitons, and some other material.
   The front and back are used at times, using pointers rpf, rpb.
   rsp is remaining space.
   powptr and comptr point to power and commutator relations.
   nng is no of new generators introduced so far.
   onng is the value of nng at the end of stage 1 (see below), which is the
   dimension of Hom-P(FM,M).
   norm is set 1 when act=1 and g is in N(P).
   separately. ipm similarly for inf4.
   stage=0,1 or 2 for H^2 and 3 or 4 for H^1. For H^2, stage=0 when computing
   Frattini module (FM), =1 when computing Hom-P(FM,M) and =2 thereafter.
   For H^1 stage=3 initially and 4 when computation of H^1 begins.
   In general tails=1 when new gens are being introduced (as tails), and the
   strings for relns in the new gens are stored at the back of rel.
   When act is true, intg and cintg are used to point to the expressions for
   the generators of the Sylow intersection Q and their conjugates, resp., as
   words in the gens of P. These are stored at the back of rel.
   subno and extno point to strings which sre the generators of the groups
   M2a and M2b, as described in Info.5. These are also stored at the back of
   rel
*/

int nqprog(void)
{
  int  i, c, **p, **q, *r, ct, oexp;
  char adn;
  if (cfm)
  /* Calculate Frattini module only  */
  {
    tails = 1;
    stage = 0;
    strcpy(outf1, outf0);
    ip = fopen(inf1, "r");
    fscanf(ip, "%d", &prime);
    fclose(ip);
    setpinv();
    calcfm(0);
    return (0);
  }
  if (rdmats() == -1)
    return (-1);
  /* If act, then we will first compute H^i(Q,M), using same algorithm as for
     H^i(P,M). First we remember the amount of space and ptrspace we had
     before we started.
  */
  if (act) {
    orsp = rsp;
    optrsp = ptrsp;
    adn = 0;
    ip = fopen(inf3, "r");
    if (ip == 0) {
      if (gap == 0)
        printf("Warning. File %s is not present.\n", inf3);
    }
    else {
      ipm = fopen(inf4, "r");
      if (ipm == 0) {
        fprintf(stderr, "Cannot open file %s.\n", inf4);
        return (-1);
      }
      while ((i = getc(ipm)) != '\n')
        ;
      inf3offset = ftell(ip);
      inf4offset = ftell(ipm);
      fclose(ip);
      fclose(ipm);
    }
  }
  while (1)
  /* If act=0 we do this only once. If act=1, each iteration in this loop
     corresponds to the action of one element g.
  */
  {
    norm = 0;
    if (act) {
      oexp = exp;
      if (ip == 0)
        i = 1;
      else if ((i = intmats()) == -1)
        return (-1);
      if (i == 1) {
        if (adn == 0)
        /* if adn=0 at this stage, there are no elements acting, so full group
           must be input, in case crel is true
        */
        {
          strcpy(inf1, inf0);
          stage = 2;
          if (ingp(1) == -1)
            return (-1);
          onng = nng;
        }
        else
          exp = oexp;
        break;
      }
      adn = 1;
    }
    /* Again we remember amount of space. If act is true, we will be storing
       some data in the back of the array rel, which will no longer be needed
       after H^i(Q,M) has been computed. rspk is recalled at end of subroutine
       spact.
    */
    if (norm) {
      inf3offset = ftell(ip);
      fclose(ip);
    }
    rspk = rsp;
    if (ch1 == 0) {
      if (norm) {
        ip = fopen(inf0, "r");
        if (ip == 0) {
          fprintf(stderr, "Cannot open file %s.\n", inf0);
          return (-1);
        }
        fscanf(ip, "%d%d%d", &prime, &exp, &facexp);
        for (i = 1; i <= 2; i++)
          while ((c = getc(ip)) != '\n')
            ;
        exp1 = -1;
        ct = 0;
        while (ct <= facexp) {
          exp1++;
          if (exp1 == exp)
            break;
          fscanf(ip, "%d", &ct);
        }
        fclose(ip);
      }
      else {
        tails = 1;
        stage = 0;
        strcpy(outf1, outf0);
        if (calcfm(matcl) == -1)
          return (-1);
        exp += nng;
      }
      if (act)
      /* Rearrange pointers intg, cintg in ptrsp, now we know how many new
         ones need to be computed.
      */
      {
        p = pcptr + ptrsp - exp1 - 1;
        q = p - exp1;
        ptrsp -= (2 * exp1);
        for (i = 1; i <= facexp; i++) {
          p[i] = intg[i];
          q[i] = cintg[i];
        }
        intg = p;
        cintg = q;
      }
    }
    if (norm) {
      tails = 0;
      stage = ch1 ? 4 : 2;
    }
    else {
      tails = 0;
      stage = ch1 ? 3 : 1;
      if (spact() == -1)
        return (-1);
      if (act == 0) {
        if (ch1 == 0) {
          subno = extno - onng;
          for (i = 1; i <= onng; i++)
            subno[i] = 0;
        }
        break;
      }
    }
    i = intact();
    if (i == -1)
      return (-1);
    if (i == 2)
      break;
    /* This ends the computation for the action of the current g. The next few
       lines prepare for the next g (if any).
    */
    rsp = orsp;
    ptrsp = optrsp;
    ip = fopen(outf2, "r");
    if (ip == 0) {
      fprintf(stderr, "Cannot open %s.\n");
      return (-1);
    }
    fscanf(ip, "%d%d%d", &prime, &dim, &ngens);
    for (i = 1; i <= ngens; i++)
      readmat(mat[i]);
    fclose(ip);
  }
  if (gap) {
    op = fopen(outfd, "w");
    if (ch1)
      fprintf(op, "COHOMOLO.CoDim1 := %d;\n", nng);
    else
      fprintf(op, "COHOMOLO.CoDim2 := %d;\n", chpdim - chsdim);
    fclose(op);
    printf("Present dimension of cohomology group is  ");
    if (ch1)
      printf("%d\n", nng);
    else
      printf("%d\n", chpdim - chsdim);
  }
  else {
    fprintf(stderr, "Present dimension of cohomology group is  ");
    if (ch1)
      fprintf(stderr, "%d\n", nng);
    else
      fprintf(stderr, "%d\n", chpdim - chsdim);
  }
  if ((ch1 && nng == 0) || (ch1 == 0 && chpdim == chsdim))
    return (2);
  if (crel)
    if (comprels() == -1)
      return (-1);
  return (0);
}

int rdmats(void)
/* reads matrices of generators of P  and set up matrix pointers.*/
{
  int i;
  int quot;
  if (act)
    ip = fopen(outf2, "r");
  else
    ip = fopen(inf2, "r");
  if (ip == 0) {
    fprintf(stderr, "Cannot open ");
    if (act)
      fprintf(stderr, "%s\n", outf2);
    else
      fprintf(stderr, "%s\n", inf2);
    return (-1);
  }
  fscanf(ip, "%d%d%d", &prime, &dim, &ngens);
  if (act == 0)
    fscanf(ip, "%d", &matcl);
  setpinv();
  quot = msp / dim;
  if (quot > mv)
    quot = mv;
  maxv = quot;
  for (i = 0; i < maxv; i++)
    vec[i] = mspace - 1 + i * dim;
  maxm = maxv / dim;
  if (maxm >= mm)
    maxm = mm - 1;
  for (i = 0; i <= maxm; i++)
    mat[i] = vec - 1 + i * dim;
  spm = mat[0];
  spv = spm[1];
  if (maxm < ngens) {
    fprintf(stderr, "Not enough mat space. Increase MSP (of MV or MM).\n");
    return (-1);
  }
  if (act == 0) {
    for (i = 1; i <= dim; i++)
      fscanf(ip, "%d", swt + i);
    for (i = 1; i <= dim; i++)
      fscanf(ip, "%d", sd1 + i);
    for (i = 1; i <= dim; i++)
      fscanf(ip, "%d", sd2 + i);
  }
  for (i = 1; i <= ngens; i++)
    readmat(mat[i]);
  fclose(ip);
  return (0);
}

int intmats(void)
/* When act=1, reads in gens for the action, computes their matrices,
   and makes a base change for this action if necesary.
*/
{
  int i, j, k, *p, l, ct, **swop, cbno2, **cbm, **newmat;
  ip = fopen(inf3, "r");
  fseek(ip, inf3offset, 0);
retry:
  fscanf(ip, "%d", &exp);
  /* if exp=0, the corresponding dcrep matrix must be skipped */
  if (exp == 0) {
    ipm = fopen(inf4, "r");
    fseek(ipm, inf4offset, 0);
    for (i = 1; i <= dim * dim; i++)
      fscanf(ipm, "%d", &j);
    inf4offset = ftell(ipm);
    fclose(ipm);
    goto retry;
  }
  if (exp == -1) {
    fclose(ip);
    return (1);
  }
  norm = exp == ngens;
  facexp = exp;
  no = exp - 1;
  for (i = 1; i <= exp; i++)
    fscanf(ip, "%d", wt + i);
  class = 1;
  for (i = 1; i <= exp; i++)
    if (wt[i] > class)
      class = wt[i];
  ct = ch1 ? exp + dim + mng : exp;
  cintg = pcptr + ptrsp - 1 - ct;
  intg = cintg - ct;
  ptrsp -= (2 * ct);
  rpb = rel + rsp - 1;
  ct = 0;
  if (norm) {
    for (i = 1; i <= exp; i++)
      if (wt[i] == 1) {
        rpb -= 3;
        intg[i] = rpb + 1;
        *(rpb + 1) = 2;
        *(rpb + 2) = i;
        *(rpb + 3) = 1;
      }
  }
  else
    for (i = 1; i <= exp; i++) {
      if (wt[i] == 1) {
        ct++;
        fscanf(ip, "%d", &l);
        p = rpb - l;
        intg[i] = p;
        *p = l;
        *cp = 0;
        while (++p <= rpb) {
          fscanf(ip, "%d", p);
          k = *p;
          p++;
          fscanf(ip, "%d", p);
          for (j = 1; j <= *p; j++)
            cp[++(*cp)] = k;
        }
        rpb -= (l + 1);
        prod(cp, mat[ct + ngens]);
      }
      while (getc(ip) != '\n')
        ;
    }
  for (i = 1; i <= exp; i++) {
    if (wt[i] == 1) {
      fscanf(ip, "%d", &l);
      p = rpb - l;
      cintg[i] = p;
      *p = l;
      while (++p <= rpb)
        fscanf(ip, "%d", p);
      rpb -= (l + 1);
    }
    while (getc(ip) != '\n')
      ;
  }
  rsp = rpb - rel + 1;
  if (norm)
    return (0);
  if (maxm < ngens + ct || maxm < 2 * exp + 3) {
    printf("Not enough mat space. Increase MSP (of MV or MM).\n");
    return (-1);
  }
  for (i = 1; i <= ct; i++) {
    swop = mat[i];
    mat[i] = mat[i + ngens];
    mat[i + ngens] = swop;
  }
  ngens = ct;
  cbno = 2 * exp + 1;
  cbm = mat[cbno];
  /* mat[cbno] is the action base change matrix */
  if (cbdef(1, ct, cbno, sd1, sd2, swt, &matcl))
    printf("No action base change.\n");
  else {
    printf("Action base change matrix:\n");
    for (i = 1; i <= dim; i++) {
      for (j = 1; j <= dim; j++)
        printf("%3d", cbm[i][j]);
      printf("\n");
    }
    cbno2 = cbno + 2;
    inv(cbm, mat[cbno + 1]);
    newmat = mat[cbno2];
    *cp = 3;
    cp[1] = cbno;
    cp[3] = cbno + 1;
    for (i = 1; i <= ct; i++) {
      cp[2] = i;
      prod(cp, newmat);
      mat[cbno2] = mat[i];
      mat[i] = newmat;
      newmat = mat[cbno2];
    }
  }
  return (0);
}

int calcfm(int steps)
/* Computes the Frattini extension to depth steps. steps=0 computes complete
   Frattini extension, but this is only used for testing.
   At each stage, the NQA is applied to the group computed so far.
*/
{
  int  i, j, k, l, d, cl, st, dp, bd, ed;
  char inp;
  printf("Computing Frattini module to depth %d.\n", steps);
  st = 1;
  depth = -1;
  inp = act ? 0 : 1;
  /* if inp=0, initial input is from inf3 (which is already open as ip);
     After first step, input is always from inf1
  */
  while (depth != steps) {
    if (ingp(inp) == -1) {
      fprintf(stderr, "Input error.\n");
      return (-1);
    }
    if (steps > 0 && depth >= steps)
      break;
    depth++;
    bd = exp + 1;
    for (dp = depth - 1; dp >= 0; dp--) {
      ed = bd - 1;
      if (dp != 0) {
        i = ed - 1;
        while (dpth[i] == dp)
          i--;
        bd = i + 1;
      }
      if (class + 1 >= mcl) {
        fprintf(stderr, "class too big. Increase MCL.\n");
        return (-1);
      }
      for (cl = class + 1; cl >= 2; cl--) {
        for (i = 1; i <= facexp; i++) {
          if (dp == 0)
            bd = i + 1;
          for (j = bd; j <= ed; j++)
            if (wt[i] == 1 && wt[j] == (cl - 1))
              if (intgen(j, i) == -1)
                return (-1);
        }
        if (dp == 0)
          for (i = 1; i <= facexp; i++)
            if (wt[i] == cl - 1)
              if (intgen(i, i) == -1)
                return (-1);
        for (i = 2; i <= facexp; i++)
          if (wt[i] > 1) {
            if (dp == 0)
              bd = i + 1;
            for (j = bd; j <= ed; j++)
              if (wt[i] + wt[j] == cl) {
                k = d1[i];
                l = d2[i];
                if (assoc(j, k, l))
                  if (subrel(j, i) == -1)
                    return (-1);
              }
          }
        for (i = 1; i <= facexp; i++)
          for (j = i + 1; j <= facexp; j++) {
            if (dp == 0)
              bd = j + 1;
            for (k = bd; k <= ed; k++)
              if (wt[i] + wt[j] + wt[k] == cl && assoc(k, j, i)) {
                if ((l = prnrel()) == 0)
                  goto nextcl;
                if (l == -1)
                  return (-1);
              }
          }
        for (i = 1; i <= facexp; i++) {
          if (dp == 0)
            bd = i;
          for (j = bd; j <= ed; j++)
            if (wt[i] + wt[j] + 1 == cl) {
              if (assoc(j, i, i)) {
                if ((l = prnrel()) == 0)
                  goto nextcl;
                if (l == -1)
                  return (-1);
              }
              if (j != i && j <= facexp && assoc(j, j, i)) {
                if ((l = prnrel()) == 0)
                  goto nextcl;
                if (l == -1)
                  return (-1);
              }
            }
        }
      nextcl:;
      }
      if (nng == 0)
        break;
    }
    if (nng == 0) {
      printf("Frattini Module Complete.\n");
      printf("Final order at depth %d was:  %d ^ %d.\n", depth - 1, prime,
             exp);
      fflush(stdout);
      break;
    }
    outgp();
    printf("Order of group at depth %d is:  %d ^ %d\n", depth, prime,
           exp + nng);
    printf("Wasted space=%d.\n\n", wsp);
    fflush(stdout);
    if (st == 1) {
      strcpy(inf1, outf1);
      if (act) {
        inp = 1;
        exp1 = exp + nng;
        inf3offset = ftell(ip);
        fclose(ip);
      }
    }
    st++;
  }
  printf("Space used, used ptrspace=%d,%d.\n", rsp - (rpb - rpf),
         pcb - pcptr);
  fflush(stdout);
  return (0);
}
