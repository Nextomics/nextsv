# include <stdio.h>
# include <stdlib.h>
# include <zlib.h>
# include <string.h>
# include <ctype.h>

# include <sys/types.h>
# include <sys/stat.h>

# include "kseq.h"
# include "tk.h"


KSEQ_INIT(gzFile, gzread) 

static int qc_usage()
{
    printf(
"\n"
"Program: longreadqc (Quality control tool for long read sequencing)\n\n"
"Usage:   longreadqc fq [options]  \n"
"Options:\n"
"    -h, --help              output this usage information\n"
"    -i, --input_file        path to the input file, which can be fasta, fastq or bam. Use this argument if only you have only one input file\n"
"    -l, --input_list_file   a file that contains the paths of all the input files. Usage this argument if you have multiple input files\n"
"    -p, --out_prefix     prefix of output files. (required)\n"
"    -d, --out_dir        output directory. (default: ./longreadqc_out/)\n"
"\n"
"For example,\n"
"longreadqc fq -i FileName.fastq -d ./FileName_longreadqc_out/ -p sample01 \n"
"longreadqc fq -l SampleName.fastqs.list -d ./SampleName_longreadqc_out/ -p sample02 \n"
);

    return 0;

}


static int qc1fastq(char * input_file, int * read_length_count, int * ptotal_num_reads, int64_t *total_num_bases, int64_t *ptotal_a_cnt, int64_t * ptotal_c_cnt, int64_t *ptotal_g_cnt, int64_t *ptotal_t_cnt, double * gc_content_count)
{
    gzFile input_fp;
    kseq_t *seq;  
    char * read_seq;
    int read_len;
    int l; 
    int i;
    double read_gc_cnt;

    input_fp = gzopen(input_file, "r");
    if(! input_fp) {
        fprintf(stderr, "ERROR! Failed to gzopen file for reading: %s", input_file);
        exit(1);
    }
    seq = kseq_init(input_fp);
    while ((l = kseq_read(seq)) >= 0)
    {
        read_seq = seq->seq.s;
        read_len = strlen(read_seq);
        *ptotal_num_reads += 1; 
        *total_num_bases += read_len;
        if (read_len <= MAX_READ_LENGTH){ 
            read_length_count[read_len] += 1; 
        }else{
            read_length_count[MAX_READ_LENGTH] += 1; 
        }
        read_gc_cnt = 0;
        for(i = 0; i < read_len; i++)
        {
            if (read_seq[i] == 'A' || read_seq[i] == 'a'){
                *ptotal_a_cnt += 1;
            }else if (read_seq[i] == 'G' || read_seq[i] == 'g' ){
                *ptotal_g_cnt += 1;
                read_gc_cnt += 1 ;
            }else if (read_seq[i] == 'C' || read_seq[i] == 'c' ){
                *ptotal_c_cnt += 1;
                read_gc_cnt += 1;
            }else if (read_seq[i] == 'T' || read_seq[i] == 't' ){
                *ptotal_t_cnt +=1;
            }
        }
        gc_content_count[ (int)(100.0 * read_gc_cnt / (double) read_len + 0.5) ] ++;
    }

    kseq_destroy(seq); 

    gzclose(input_fp);

    return 0;
}

static int qc_fastqs (STRING_LIST * input_file_list, const char * full_out_prefix)
{
    int i;
    char * input_file = NULL;
    int total_num_reads = 0;
    int64_t total_num_bases = 0;
    int * read_length_count;
    int n50_read_length = 0;
    int n95_read_length = 0;
    int mean_read_length = 0;
    int median_read_length = 0;
    int n05_read_length = 0;
    int longest_read_length = 0;


    double gc_content_count[101] = {0};

    char * basic_info_file;
    FILE * basic_info_fp;
    char * read_length_histo_file1;
    char * read_length_histo_file2;
    char * read_length_histo_file3;
    char * gc_content_file;

    read_length_count = (int *) calloc (MAX_READ_LENGTH+1,sizeof(int));

    basic_info_file = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));
    read_length_histo_file1 = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));
    read_length_histo_file2 = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));
    read_length_histo_file3 = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));
    gc_content_file = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));

    sprintf(basic_info_file, "%s_basic_info.txt", full_out_prefix);
    sprintf(read_length_histo_file1, "%s_read_length_histo_binsize100.txt", full_out_prefix);
    sprintf(read_length_histo_file2, "%s_read_length_histo_binsize1k.txt", full_out_prefix);
    sprintf(read_length_histo_file3, "%s_read_length_histo_binsize10k.txt", full_out_prefix);
    sprintf(gc_content_file, "%s_gc_content.txt", full_out_prefix);

    basic_info_fp = fopen(basic_info_file, "w");

    int64_t total_a_cnt = 0; 
    int64_t total_c_cnt = 0; 
    int64_t total_g_cnt = 0; 
    int64_t total_t_cnt = 0; 
    for (i = 0; i < input_file_list->size; i++) {
        input_file = input_file_list->data_list[i];
        qc1fastq(input_file, read_length_count, &total_num_reads, &total_num_bases, &total_a_cnt, &total_c_cnt, &total_g_cnt, &total_t_cnt, gc_content_count);
    }

    get_read_length_statistics(read_length_count, MAX_READ_LENGTH, &n50_read_length, &mean_read_length, &median_read_length, &n05_read_length, &n95_read_length, &longest_read_length);

    fprintf (basic_info_fp, "output prefix\t%s\n", full_out_prefix);
    fprintf (basic_info_fp, "total number of reads\t%d\n", total_num_reads);
    fprintf (basic_info_fp, "total number of bases\t%lld\n", total_num_bases);

    fprintf (basic_info_fp, "longest read length\t%d\n", longest_read_length);
    fprintf (basic_info_fp, "N50 read length\t%d\n", n50_read_length);
    fprintf (basic_info_fp, "mean read length\t%d\n", mean_read_length);
    fprintf (basic_info_fp, "median read length\t%d\n", median_read_length);
    fprintf (basic_info_fp, "N05 read length\t%d\n", n05_read_length);
    fprintf (basic_info_fp, "N95 read length\t%d\n", n95_read_length);
    fprintf (basic_info_fp, "%%GC\t%.2f\n", ((double) (total_g_cnt + total_c_cnt) / (double) total_num_bases * 100.0));

    // basic statistics
    // gc content 
    // read length histogram, cumulative read length, 
    // base quality distribution
    output_read_length_histo (read_length_count, MAX_READ_LENGTH, read_length_histo_file1, 100);
    output_read_length_histo (read_length_count, MAX_READ_LENGTH, read_length_histo_file2, 1000);
    output_read_length_histo (read_length_count, MAX_READ_LENGTH, read_length_histo_file3, 10000);
    //output_gc_content_histo (read_length_count, MAX_READ_LENGTH, read_length_histo_file1);
    for (i = 0; i <= 100; i++) {
        //fprintf(stdout, "%d\t%.f\n", i, gc_content_count[i]);
    }

    fclose(basic_info_fp);

    return 0;
}

int main_qcfq(int argc, char * argv[])
{
    char * input_file = NULL; 
    char * input_list_file  = NULL;
    char * out_prefix = NULL;
    char * full_out_prefix = NULL;
    char * line = NULL;
    char * out_dir = NULL;
    STRING_LIST * input_file_list;
    line = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));
    out_dir = (char *) calloc(4096, sizeof(char)); 

    if (argc < 2) { 
        qc_usage(); 
        return 1; 
    }

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 ) {
        qc_usage();
        return 0;
    }

    input_file_list = init_string_list(1000, MAX_PATH_LENGTH);

    int index = 1;
    while (index < argc){

        if (strcmp(argv[index], "--input_file") == 0 || strcmp(argv[index], "-i" ) == 0){
            if (index +1 >= argc){ qc_usage(); return 1;};
            input_file = argv[index+1];
            index += 2;
        }else if (strcmp(argv[index], "--input_list_file") == 0 || strcmp(argv[index], "-l" ) == 0){
            if (index +1 >= argc){ qc_usage(); return 1;};
            input_list_file = argv[index+1];
            index += 2;
        }else if (strcmp(argv[index], "--out_dir") == 0 || strcmp(argv[index], "-d" ) == 0){
            if (index +1 >= argc){ qc_usage(); return 1;};
            sprintf(out_dir, "%s", argv[index+1]);
            index += 2;
        }else if (strcmp(argv[index], "--out_prefix") == 0 || strcmp(argv[index], "-p" ) == 0){
            if (index +1 >= argc){ qc_usage(); return 1;};
            out_prefix = argv[index+1];
            index += 2;
        }else{
            qc_usage(); 
            return 1;
        }
    }

    if (input_file != NULL && input_list_file != NULL){
        fprintf(stderr, "--input_file and --input_list_file cannot be specified together");
        qc_usage();
        exit(1);
    }
    if (input_file == NULL && input_list_file == NULL){
        fprintf(stderr, "ERROR! No input files. Both --input_file and --input_list_file were not specified.");
        qc_usage();
        exit(1);
    }
    if (out_dir[0] == 0) { 
        sprintf(out_dir, "./longreadqc_out/");
    }
    if (out_prefix == NULL) {
        fprintf(stderr, "ERROR! --out_prefix was not specified.");
        qc_usage();
        exit(1);

    }
    struct stat st = {0};
    if (stat(out_dir, &st) == -1) {
        mkdir(out_dir, 0755);
    }

    if (input_file != NULL){
        add1input_file(input_file_list, input_file);
    }
    if (input_list_file != NULL){
        FILE * input_list_fp;
        input_list_fp = fopen(input_list_file, "r");
        if(input_list_fp == NULL){
            fprintf(stderr, "ERROR! Failed to open file for reading: %s", input_list_file);
            exit(1);
        }
        while (fgets(line, MAX_PATH_LENGTH, input_list_fp) != NULL)
        {
            add1input_file(input_file_list, line);
        }

    }

    int i;
    fprintf(stdout, "%d input file(s):\n", input_file_list->size);
    for (i = 0; i < input_file_list->size; i++)
    {
        fprintf(stdout, "%s\n", input_file_list->data_list[i]);
    }
    fprintf(stdout, "output directory: %s\n", out_dir);

    full_out_prefix = join_path(out_dir, out_prefix);

    qc_fastqs(input_file_list, full_out_prefix);

    free(out_dir);

    return 0;

}
