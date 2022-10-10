#include "permfns.h"
#include "defs.h"

extern short npt, cp[], orb[], *pptr[], pno[];
extern FILE *ip, *op;

short orbitsv(short pt, short * sv, short lo)
/* Computes orbit of pt under perms listed in pno, and writes
   Schreier vector into sv
*/
{
  short u, v, w, x, y, z;
  if (lo == 0) {
    for (u = 1; u <= npt; u++)
      sv[u] = 0;
    orb[1] = pt;
    lo = 1;
    sv[pt] = -1;
  }
  for (x = 1; x <= lo; x++) {
    z = orb[x];
    for (y = 1; y <= *pno; y++) {
      w = pno[y];
      v = pptr[w][z];
      if (sv[v] == 0) {
        lo++;
        orb[lo] = v;
        sv[v] = w + 1;
      }
    }
  }
  return (lo);
}

short addsv(short pt, short * sv)
/* The Schreier vector sv is applied to the point pt, and the resulting word
   added to the end of cp. It is assumed that pt is in the orbit. If not, the
   program will break down.
   Externals: cp.
*/
{
  short pn;
  pn = sv[pt];
  while (pn != -1) {
    (*cp)++;
    cp[*cp] = pn;
    pt = pptr[pn][pt];
    pn = sv[pt];
  }
  return (0);
}

short image(short pt)
/* The image of pt under cp is computed and returned.
   Externals: pptr,cp.
*/
{
  short i;
  for (i = 1; i <= *cp; i++)
    pt = pptr[cp[i]][pt];
  return (pt);
}

short invert(short * ptr1, short * ptr2)
/* permutation ptr1 is inverted and put in ptr2.
   Externals: npt.
*/
{
  short i;
  for (i = 1; i <= npt; i++)
    ptr2[ptr1[i]] = i;
  return (0);
}

short readperm(short * ptr)
/* The next npt numbers from input ip are read. These should form a
   permutation on 1,2,3,...,npt. If not  2  is returned. If the perm is the
   identity, 1 is returned, otherwise 0.
   Externals:  npt,orb,ip.
*/
{
  short i, j, id;
  id = 1;
  for (i = 1; i <= npt; i++)
    orb[i] = 0;
  for (i = 1; i <= npt; i++) {
    fscanf(ip, "%hd", ptr + i);
    j = ptr[i];
    if (j <= 0 || j > npt || orb[j]) {
      fprintf(stderr, "perm[%d]=%d\n", i, j);
      return (2);
    }
    orb[j] = 1;
    if (id && j != i)
      id = 0;
  }
  return (id);
}

short printvec(short * ptr, short e)
/* Points ptr[1] to ptr[npt+e] are output to op, followed by new line. The
   first npt of these will be a permutation or a Schreier vector. e=0 or 1.
  Externals: npt,op.
*/
{
  short i;
  if (npt >= 10000)
    for (i = 1; i <= npt; i++)
      fprintf(op, "%6d", ptr[i]);
  else if (npt >= 1000)
    for (i = 1; i <= npt; i++)
      fprintf(op, "%5d", ptr[i]);
  else
    for (i = 1; i <= npt; i++)
      fprintf(op, "%4d", ptr[i]);
  fprintf(op, "    ");
  for (i = npt + 1; i <= npt + e; i++)
    fprintf(op, "%4d", ptr[i]);
  fprintf(op, "\n");
  return (0);
}

short readvec(short * ptr, short e)
/* The next npt+e points from ip are read into array ptr.
   Externals: npt,ip.
*/
{
  short i;
  for (i = 1; i <= npt + e; i++)
    fscanf(ip, "%hd", ptr + i);
  return (0);
}

short readbaselo(short nb, short * base, short * lorb)
/* The nb base points are read into base, and the nb orbit lengths into lorb,
   from ip.
   Externals: ip.
*/
{
  short i;
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%hd", base + i);
  for (i = 1; i <= nb; i++)
    fscanf(ip, "%hd", lorb + i);
  return (0);
}

short printbaselo(short nb, short * base, short * lorb)
/* base and lorb are printed to op.
   Externals: op.
*/
{
  short i;
  if (npt >= 1000) {
    for (i = 1; i <= nb; i++)
      fprintf(op, "%5d", base[i]);
    fprintf(op, "\n");
    for (i = 1; i <= nb; i++)
      fprintf(op, "%5d", lorb[i]);
    fprintf(op, "\n");
  }
  else {
    for (i = 1; i <= nb; i++)
      fprintf(op, "%4d", base[i]);
    fprintf(op, "\n");
    for (i = 1; i <= nb; i++)
      fprintf(op, "%4d", lorb[i]);
    fprintf(op, "\n");
  }
  return (0);
}

short printpsv(short nb, short * gno, short ** svptr)
/* Permutationsnos gno[1],...,gno[*gno] are output (up to npt+1), and then
   Schreier vectors svptr[1],...,svptr[nb] are output to op.
   Externals: npt,pptr,orb,op.
*/
{
  short i, j, k;
  for (i = 1; i <= *gno; i++) {
    j = gno[i];
    orb[j + 1] = 2 * i - 1;
    printvec(pptr[j], 1);
  }
  for (i = 1; i <= nb; i++) {
    for (j = 1; j <= npt; j++) {
      k = svptr[i][j];
      if (k > 0)
        k = orb[k];
      fprintf(op, "%4d", k);
    }
    fprintf(op, "\n");
  }
  return (0);
}

short readpsv(short e, short nb, short nperms, short ** svptr)
/* nperms permutations (up to npt+1) are read into perm nos
   pptr[e],pptr[e+2],. and pptr[e+2x] is inverted into pptr[e+2x+1], for
   x=1,...,nperms. Then the Screier vectors svptr[1],...,svptr[nb] are read
   from ip. Externals: pptr,npt,ip.
*/
{
  short i, j, *k;
  for (i = 1; i <= nperms; i++) {
    j = e + 2 * i - 2;
    readvec(pptr[j], 1);
    invert(pptr[j], pptr[j + 1]);
  }
  for (i = 1; i <= nb; i++) {
    readvec(svptr[i], 0);
    for (j = 1; j <= npt; j++) {
      k = svptr[i] + j;
      if (*k > 0)
        *k += e;
    }
  }
  return (0);
}

void seeknln(void)
{
  while (getc(ip) != '\n')
    ;
}
