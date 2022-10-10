#include "defs.h"

extern short npt, nb, npt1, exp, prime, pinv[], cp[], power[], wt[], base[],
    ngno[], igno[], *svptr[], *pptr[], d1[];

/* Let the central series for P be P=P(n)>...>P(1)>P(0)=1. The procedures
   in this file express a permutation g in P in in its normal form in the
   PCP generators. Let the Schreier vector generators be g(n),...,g(1) and the
   PCP generators h(n),...,h(1), with g(i),h(i) in P(i)-P(i-1).
   Then we have g(i)=h(i)^x * k, with k in P(i-1), and the exponent x is
   stored as power[i].
*/
void firstgen(short * p, short * hg, short * co)
/* Let the perm p be  g(i)^t * k with k in H(i-1). This procedure computes i
   and t, and returns them as *hg and *co.
*/
{
  short i, j;
  *cp = 0;
  for (i = 1; i <= nb; i++)
    addsv(image(p[base[i]]), svptr[i]);
  *co = 0;
  *hg = 0;
  for (i = 1; i <= *cp; i++) {
    j = igno[cp[i]];
    if (j > *hg) {
      *hg = j;
      *co = 1;
    }
    else if (j == *hg)
      (*co)++;
  }
  (*co) %= prime;
}

int express(short * p, short * relc, int nwt)
/* This computes the full PCP expression for the perm p, and puts as a gen-pow
   string in rel, preceded by its length. If weights are involved, then nwt is
   the weight of p at this point. Otherwise nwt=0.
*/
{
  short hgen, h, l, coeff, pow, m, n, *ip, *gp, pt;
  l = 0;
  h = exp;
  ip = p + npt1;
  while (1) {
    firstgen(p, &hgen, &coeff);
    if (hgen > h) {
      fprintf(stderr, "Expression error.\n");
      return (-1);
    }
    if (hgen == 0)
      break;
    l++;
    relc[l] = exp + 1 - hgen;
    if ((nwt && wt[hgen] < nwt) ||
        (nwt && wt[hgen] == nwt && d1[hgen] == 0)) {
      power[hgen] = pinv[coeff];
      l++;
      relc[l] = 1;
      *relc = l;
      return (hgen);
    }
    /* We return at this point, because PCP gen no hgen will be replaced by p
     */
    pow = (coeff * power[hgen]) % prime;
    l++;
    relc[l] = pow;
    gp = pptr[ngno[hgen]];
    for (n = 1; n <= npt; n++) {
      pt = ip[n];
      for (m = 1; m <= pow; m++)
        pt = gp[pt];
      p[pt] = n;
    }
    invert(p, ip);
    h = hgen - 1;
  }
  *relc = l;
  return (0);
}

void setpinv(void)
{
  short i, j;
  for (i = 0; i < prime; i++)
    pinv[i] = 0;
  for (i = 1; i < prime; i++)
    if (pinv[i] == 0)
      for (j = 1; j < prime; j++)
        if (i * j % prime == 1) {
          pinv[i] = j;
          pinv[j] = i;
          break;
        }
}
