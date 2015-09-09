#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "fermi.h"
#include "mag.h"

//modified from Heng Li's example.c of fermi software.

int main(int argc, char *argv[])
{
	int c, ec_k = -1, unitig_k = -1;
	char *seq, *qual;
	FILE *fp;
	int64_t l;
	while ((c = getopt(argc, argv, "k:l:")) >= 0) {
		switch (c) {
			case 'k': ec_k = atoi(optarg); break;
			case 'l': unitig_k = atoi(optarg); break;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "Local assembler for small peak regions. Output cleaned unitigs.\nUsage: SAPPER_fermi [-k ecKmer] [-l utgKmer] <in.fq> <out.fq>\n");
		return 1;
	}
	l = fm6_api_readseq(argv[optind], &seq, &qual);
	fp = fopen( argv[optind+1], "w");
	fm6_api_correct(ec_k, l, seq, qual);
	mag_t *g;
	free(qual);
	g = fm6_api_unitig(unitig_k, l, seq);
	magopt_t *opt = mag_init_opt();
	//opt->flag |= MOG_F_AGGRESSIVE | MOG_F_CLEAN;
	opt->flag |= MOG_F_CLEAN;
	mag_g_clean(g, opt);
	free(opt);
	mag_g_fprint(g, fp);
	mag_g_destroy(g);
	free(seq);
	return 0;
}
