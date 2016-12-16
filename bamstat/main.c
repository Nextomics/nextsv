#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <sam.h>


int bamstat(char * bam)
{	
	samFile * in;
	bam_hdr_t * hdr;
	bam1_t * b;
	bam1_core_t * c;

	int64_t i;
	int64_t tot_base_c;
	int64_t tot_read_c;
	int64_t mapped_read_c;
	int64_t mapped_base_c;
	int64_t second_aln_c;
	int64_t suppl_aln_c;
	int ret;

    b = bam_init1(); 
	c = &(b->core);
	tot_base_c = tot_read_c = 0;
	mapped_read_c = mapped_base_c = 0;
	second_aln_c = suppl_aln_c = 0;	

	in = sam_open(bam, "rb");
	if (NULL == in){
		fprintf (stderr, "ERROR: failed to open file: %s\n", bam);
		abort();
	}
	hdr = sam_hdr_read(in);
	while (ret = sam_read1(in, hdr, b) >= 0){ 
		if (c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY) 
		{
			if (c->flag & BAM_FSECONDARY){
				second_aln_c++;
			}
			if (c->flag & BAM_FSUPPLEMENTARY){
				suppl_aln_c++;
			}
			continue;
		}

		tot_read_c++;
		tot_base_c += c->l_qseq;

		if (c->flag & BAM_FUNMAP){
			continue;
		}

		mapped_read_c++;
		mapped_base_c += bam_endpos(b) - c->pos;

	}

	printf ("Total read count (M):  \t%.6f\n", (double)tot_read_c/1e6);
	printf ("Total base count (G):  \t%.6f\n", (double)tot_base_c/1e9);
	printf ("Mapped read count (M): \t%.6f\n", (double)mapped_read_c/1e6);
	printf ("Mapped base count (M): \t%.6f\n", (double)mapped_base_c/1e6);
	printf ("\n");

	return 0;
}

int usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bamstat (mapping statistics for bam files)\n"); 
	fprintf(stderr, "Version: 0.4.0\n"); 
	fprintf(stderr, "Contact: Li Fang (fangli@grandomics.com)\n"); 
	fprintf(stderr, "Usage:   bamstat in.bam > out.stat.txt\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char * argv[])
{
	if (argc < 2) {
		return usage();
	}
	
	bamstat(argv[1]);

	return 0;
}

