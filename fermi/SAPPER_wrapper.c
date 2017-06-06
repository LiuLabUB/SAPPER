#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "fermi.h"
#include "mag.h"

//modified from Heng Li's example.c of fermi software.

void assemble( char * infq, char * outfq, int unitig_k, int merge_min_len )
{
  int ec_k = -1;
  char *seq, *qual;
  FILE *fp;
  int64_t l;
  fp = fopen( outfq, "wb");
  l = fm6_api_readseq(infq, &seq, &qual);
  fm6_api_correct(ec_k, l, seq, qual);
  mag_t *g;
  free(qual);
  g = fm6_api_unitig(unitig_k, l, seq);
  magopt_t *opt = mag_init_opt();
  opt->flag |= MOG_F_CLEAN;
  opt->min_merge_len = merge_min_len;
  mag_g_clean(g, opt);
  free(opt);
  mag_g_fprint(g, fp);
  mag_g_destroy(g);
  free(seq);
  fclose(fp);
  return;
}
