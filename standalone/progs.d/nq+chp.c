#include "defs.h"

extern char    inf1[], outf[80];
extern short **mat[], pinv[], **spm, *spv, prime, dim, cp[], exp, opmats;
FILE *         ip, *ips, *op;

int matact(int inter)
{
  short substeps, subexp, subdp, *commno, *commvec, **commer, **commkeep,
      intexp, ncommer, ncomm, i, j, k, l, m, n, *ptr, *ptre, *ptr2, incr;
  char subfile[4], infc[80];
  commer = mat[3 * exp + 1];
  commkeep = mat[3 * exp + 2];
  commno = spm[1];
  commvec = spm[2];
  if (inter) {
    if ((ip = fopen(inf1, "r")) == 0) {
      fprintf(stderr, "Cannot open %s\n", inf1);
      return (-1);
    }
    strcpy(subfile, "axa");
    incr = exp;
  }
  else {
    incr = 0;
    intexp = exp;
    strcpy(subfile, "xa");
  }

  /* The following loop is executed once only if inter = 0 (Checking relations
     for P). If inter=1, it is executed once for each Sylow intersection Q in
     the file inf1 (=gpname.sc?).
  */
  while (1) {
    if (inter) {
      fscanf(ip, "%hd", &intexp);
      if (intexp == 0) {
        subfile[0]++;
        subfile[2] = 'a';
        continue;
      }
      if (intexp == -1)
        break;
      if (intexp >= exp) {
        fprintf(stderr, "Sylow intersection too large.\n");
        return (-1);
      }
      seeknln();
      seeknln();
      /* Read generators of Q and compute matrices */
      for (i = 1; i <= intexp; i++) {
        fscanf(ip, "%hd", &l);
        l /= 2;
        ptr = ptre = cp;
        for (j = 1; j <= l; j++) {
          fscanf(ip, "%hd%hd", &m, &n);
          if (m <= 0 || m > exp || n <= 0 || n > prime) {
            fprintf(stderr, "Invalid data in %s.\n", inf1);
            return (-1);
          }
          ptre += n;
          while (++ptr <= ptre)
            *ptr = m;
          ptr = ptre;
        }
        *cp = ptre - cp;
        prod(cp, mat[exp + i]);
      }
      if (opmats) {
        strcpy(outf, inf1);
        strncat(outf, subfile, 1);
        strcat(outf, "mats");
        op = fopen(outf, "w");
        fprintf(op, "%4d%4d%4d\n", prime, dim, intexp);
        for (i = 1; i <= intexp; i++)
          printmat(mat[exp + i]);
        fclose(op);
      }
      /* Skip rest of data for this Q */
      k = 3 + intexp * (intexp + 1) / 2;
      for (i = 1; i <= k; i++)
        seeknln();
    }

  /* The following loop (restart) is executed once for each file of inf1xa,
     inf1xb,... of subgroups of Q.
  */
  restart:
    strcpy(infc, inf1);
    strcat(infc, subfile);
    if ((ips = fopen(infc, "r")) == 0)
    /* No more files of subgroups. */
    {
      if (inter == 0)
        break;
      subfile[0]++;
      subfile[2] = 'a';
      continue;
    }
    fscanf(ips, "%hd%hd", &subexp, &substeps);
    if (subexp <= 0 || subexp > intexp) {
      fprintf(stderr, "Invalid data in %s.\n", infc);
      return (-1);
    }
    printf("File %s:  substeps=%d.\n", infc, substeps);
    for (i = 1; i <= subexp; i++) {
      if (subexp == intexp) {
        *cp = 1;
        cp[1] = i + incr;
      }
      else {
        fscanf(ips, "%hd", &l);
        l /= 2;
        ptr = ptre = cp;
        for (j = 1; j <= l; j++) {
          fscanf(ips, "%hd%hd", &m, &n);
          if (m <= 0 || m > intexp || n <= 0 || n > prime) {
            fprintf(stderr, "Invalid data in %s.\n", infc);
            return (-1);
          }
          ptre += n;
          while (++ptr <= ptre)
            *ptr = m + incr;
          ptr = ptre;
        }
        *cp = ptre - cp;
      }
      prod(cp, mat[2 * exp + i]);
    }

    for (i = 1; i <= dim; i++) {
      for (j = 1; j <= dim; j++)
        commer[i][j] = 0;
      commer[i][i] = 1;
    }
    ncommer = dim;
    for (subdp = 1; subdp <= substeps; subdp++) {
      ncomm = 0;
      for (i = 1; i <= dim; i++)
        commno[i] = 0;
      for (i = 1; i <= ncommer; i++)
        for (j = 1; j <= subexp; j++) {
          k = comm(commer[i], commvec, mat[2 * exp + j]);
          if (subdp == substeps) {
            if (k != 0) {
              fprintf(stderr, "Error: subdp,i,j,  inf=%d,%d,%d,%s\n", subdp,
                      i, j, infc);
              subfile[1 + inter]++;
              fclose(ips);
              goto restart;
            }
          }
          else {
            ptr = commvec + dim + 1;
            for (k = dim; k >= 1; k--)
              if ((l = *(--ptr)) != 0) {
                n = pinv[l];
                for (m = k; m >= 1; m--) {
                  *ptr *= n;
                  *(ptr--) %= prime;
                }
                if ((n = commno[k]) != 0) {
                  ptr = commvec + k;
                  ptr2 = commkeep[commno[k]] + k;
                  while (ptr > commvec) {
                    *(ptr) -= *ptr2;
                    if (*ptr < 0)
                      *ptr += prime;
                    ptr--;
                    ptr2--;
                  }
                  ptr = commvec + k;
                }
                else {
                  ncomm++;
                  commno[k] = ncomm;
                  ptr = commkeep[ncomm];
                  ptre = ptr + dim;
                  ptr2 = commvec;
                  while (++ptr <= ptre)
                    *ptr = *(++ptr2);
                  break;
                }
              }
          }
        }

      if (subdp < substeps) {
        fscanf(ips, "%hd", &i);
        if (i != 0) {
          subexp = i;
          if (subexp <= 0 || subexp > intexp) {
            fprintf(stderr, "Invalid data in %s.\n", infc);
            return (-1);
          }
          for (i = 1; i <= subexp; i++) {
            if (subexp == intexp) {
              *cp = 1;
              cp[1] = i + incr;
            }
            else {
              fscanf(ips, "%hd", &l);
              l /= 2;
              ptr = ptre = cp;
              for (j = 1; j <= l; j++) {
                fscanf(ips, "%hd%hd", &m, &n);
                if (m <= 0 || m > intexp || n <= 0 || n > prime) {
                  fprintf(stderr, "Invalid data in %s.\n", infc);
                  return (-1);
                }
                ptre += n;
                while (++ptr <= ptre)
                  *ptr = m + incr;
                ptr = ptre;
              }
              *cp = ptre - cp;
            }
            prod(cp, mat[2 * exp + i]);
          }
        }
        for (i = 1; i <= ncomm; i++) {
          ptr = commer[i];
          ptre = ptr + dim;
          ptr2 = commkeep[i];
          while (++ptr <= ptre)
            *ptr = *(++ptr2);
        }
        ncommer = ncomm;
      }
    }
    printf("File %s is ok.\n", infc);
    fclose(ips);
    subfile[1 + inter]++;
    goto restart;
  }
}

void seeknln(void)
{
  while (getc(ip) != '\n')
    ;
}
