# include <stdio.h>
# include <stdlib.h>
# include <zlib.h>
# include <string.h>
# include <ctype.h>
# include "tk.h"
# include "kseq.h"

KSEQ_INIT(gzFile, gzread) 

static int qc_usage()
{
    printf(
"\n"
"Program: longreadqc (Quality control tool for long read sequencing)\n\n"
"Usage:   longreadqc fq \n"
"Options:\n"
"    -h, --help              output this usage information\n"
"    -i, --input_file        path to the input file, which can be fasta, fastq or bam. Use this argument if only you have only one input file\n"
"    -l, --input_list_file   a file that contains the paths of all the input files. Usage this argument if you have multiple input files\n"
"    -t, --input_type        format of the input file(s), (`fasta` or `fastq` or `bam`). (infered from filename extension)\n"
"    -o, --out_prefix        prefix of the output files. (default: ./InputFileName.longreadqc or ./InputListFileName.longreadqc)\n");

    return 0;

}


int main_qcfa(int argc, char * argv[])
{

    return 0;

}
