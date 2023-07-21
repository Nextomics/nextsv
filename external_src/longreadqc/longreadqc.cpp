# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "tk.h"

int main_qcfq  (int argc, char * argv[]);
int main_qcfa  (int argc, char * argv[]);
int main_qcbam (int argc, char * argv[]);
int main_qcpaf (int argc, char * argv[]);
int main_filterfq (int argc, char * argv[]);

int usage()
{
    printf(
"\n"
"Program: longreadqc (Quality control tool for long read sequencing)\n"
"\n"
"Usage:   longreadqc fq          quality control for FASTQ files\n"
"         longreadqc filterfq    filter FASTQ reads \n"
"         longreadqc paf         quality control for PAF file \n"
//"         longreadqc fa          quality control for fasta files\n"
//"         longreadqc bam         quality control for bam files\n"
);
    return 0;
}

int main(int argc, char * argv[])
{
    if (argc < 2) { 
        usage(); 
        return 1; 
    }

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0) {
        usage();
        return 0;
    }

    int ret = 0;

    if(strcmp(argv[1], "fq") == 0){
        ret = main_qcfq(argc-1, argv+1);
    }else if (strcmp(argv[1], "fa") == 0){
        ret = main_qcfa(argc-1, argv+1);
    }else if (strcmp(argv[1], "paf") == 0){
        ret = main_qcpaf(argc-1, argv+1);
    }else if (strcmp(argv[1], "bam") == 0){
        ret = main_qcbam(argc-1, argv+1);
    }else if (strcmp(argv[1], "filterfq") == 0){
        ret = main_filterfq(argc-1, argv+1);
    }else{
        usage();
        return 0;
    }
    return ret;
}
