#include "defs.h"
#include "permfns.h"

extern short npt, pno[], expcp[], cp[], genorb[], *pptr[], *expptr[];

int expandp(int nb, short * p, short * base, short ** svptr)
/* Generates the Schreier vector expression for perm p into cp.
   Externals: cp.
*/
{
  short i, j, *k;
  *cp = 0;
  for (i = 1; i <= nb; i++) {
    j = image(p[base[i]]);
    k = svptr[i];
    if (k[j] == 0)
      return (2);
    addsv(j, k);
  }
  return (0);
}

int resetsv(int nb, short * base, short * lorb, short * gno, short ** svptr)
/* Recomputes the orbits and Schreier vectors of perm nos
   gno[1],...,gno[*gno]. Externals: pno,pptr,npt.
*/
{
  short i, j, npt1;
  npt1 = npt + 1;
  for (i = 1; i <= nb; i++) {
    *pno = 0;
    for (j = 1; j <= *gno; j++)
      if (pptr[gno[j]][npt1] >= i) {
        (*pno)++;
        pno[*pno] = gno[j];
        lorb[i] = orbitsv(base[i], svptr[i], 0);
      }
  }
  return (0);
}

int allorbs(short * lorb, short * orno)
/* Computes all orbits of perm nos pno[1],...,pno[*pno]. Point i lies in orbit
   no orno[i], which has length lorb[orno[i]]. lorb[0)=no of orbits. The
   extern array genorb is used to list points in current orbit.
   Externals: genorb,npt,pno,pptr.
*/
{
  short orct, lo, u, v, w, x, y, z;
  for (u = 1; u <= npt; u++)
    orno[u] = 0;
  orct = 0;
  for (u = 1; u <= npt; u++)
    if (orno[u] == 0) {
      orct++;
      orno[u] = orct;
      lo = 1;
      genorb[1] = u;
      for (x = 1; x <= lo; x++) {
        z = genorb[x];
        for (y = 1; y <= *pno; y++) {
          w = pno[y];
          v = pptr[w][z];
          if (orno[v] == 0) {
            lo++;
            genorb[lo] = v;
            orno[v] = orct;
          }
        }
      }
      lorb[orct] = lo;
    }
  *lorb = orct;
  return (0);
}

int backimage(int pt)
/* This computes and returns the image of pt under the inverse of the perm cp.
   Externals: cp,pptr.
*/
{
  short i;
  for (i = *cp; i >= 1; i--)
    pt = pptr[cp[i] - 1][pt];
  return (pt);
}

int exprep(int pt, int no, short * sv)
/* The word for pt is computed using Schreier vector sv, and the corresponding
   perm (i.e. the inverse of cp) stored in expptr[no].
   Externals: npt,expptr,cp.
*/
{
  short i;
  *cp = 0;
  addsv(pt, sv);
  for (i = 1; i <= npt; i++)
    expptr[no][i] = backimage(i);
  return (0);
}

int expimage(int pt)
/* The image of pt is computed and returned under the word stored in expcp.
   Externals: expcp,expptr.
*/
{
  short i;
  for (i = *expcp; i >= 1; i--)
    pt = expptr[expcp[i]][pt];
  return (pt);
}
