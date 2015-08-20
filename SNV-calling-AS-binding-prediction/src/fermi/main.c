#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include "fermi.h"

int main_splitfa(int argc, char *argv[]);
int main_fltuniq(int argc, char *argv[]);
int main_cg2cofq(int argc, char *argv[]);
int main_pe2cofq(int argc, char *argv[]);
int main_trimseq(int argc, char *argv[]);

int main_chkbwt(int argc, char *argv[]);
int main_unpack(int argc, char *argv[]);
int main_exact(int argc, char *argv[]);
int main_merge(int argc, char *argv[]);
int main_sub(int argc, char *argv[]);
int main_build(int argc, char *argv[]);
int main_correct(int argc, char *argv[]);
int main_unitig(int argc, char *argv[]);
int main_clean(int argc, char *argv[]);
int main_cnt2qual(int argc, char *argv[]);
int main_seqsort(int argc, char *argv[]);
int main_remap(int argc, char *argv[]);
int main_scaf(int argc, char *argv[]);
int main_contrast(int argc, char *argv[]);
int main_bitand(int argc, char *argv[]);
int main_recode(int argc, char *argv[]);

int main_ropebwt(int argc, char *argv[]);
int main_example(int argc, char *argv[]);

double rssmem();
double cputime();
double realtime();
void liftrlimit();
/*
#include "rld.h"
int main_test(int argc, char *argv[])
{
	extern uint64_t fm_multi_backward_search(int n, rld_t *const*e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end);
	int i, l;
	uint64_t beg, end;
	l = strlen(argv[1]);
	for (i = 0; i < l; ++i)
		argv[1][i] = seq_nt6_table[(int)argv[1][i]];
	if (argc > 3) {
		int n = argc - 2;
		rld_t **e;
		e = alloca(n * sizeof(void*));
		for (i = 2; i < argc; ++i)
			e[i - 2] = rld_restore(argv[i]);
		fm_multi_backward_search(n, e, l, (uint8_t*)argv[1], &beg, &end);
		printf("%lld,%lld\n", beg, end);
	} else {
		rld_t *e;
		e = rld_restore(argv[2]);
		fm_backward_search(e, l, (uint8_t*)argv[1], &beg, &end);
		printf("%lld,%lld\n", beg, end);
	}
	return 0;
}
*/
int main(int argc, char *argv[])
{
	//int i, ret = 0;
	int ret = 0;
	//double start = realtime();
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: SNVAS_fermi\n");
		fprintf(stderr, "Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n\n");
		fprintf(stderr, "Command: SNVAS_fermi [-ceU] [-k ecKmer] [-l utgKmer] <in.fq>\n\n");
		fprintf(stderr, "\n");
		return 1;
	}

	ret = main_example(argc, argv);

	/*if (ret == 0 && fm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, FERMI_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; RSS: %.3f MB\n", __func__, realtime() - start, cputime(), rssmem());
	}*/
	return ret;
}
