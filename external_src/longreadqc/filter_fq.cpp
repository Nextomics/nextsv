# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <zlib.h>
# include <string.h>
# include <ctype.h>

# include <sys/types.h>
# include <sys/stat.h>
# include <iostream>
# include <string>
# include <unordered_map>

# include "tk.h"


char is_letter[256] = {0};
char is_qual_value[256] = {0};

static int filterfq_usage()
{
    printf(
"\n"
"Program: longreadqc (Quality control tool for long read sequencing)\n\n"
"Usage:   longreadqc filterfq [options]  \n"
"Options:\n"
"    -h, --help              output this usage information\n"
"    -i, --input_file        path to the input fastq file. Use this argument if only you have only one input file\n"
"    -l, --input_list_file   a file that contains the paths of all the input fastq files. Usage this argument if you have multiple input files\n"
"    -p, --out_prefix        prefix of the output file. (required)\n"
"    -n, --num_split_file    split the output file to n files. (not required, default: n=1)\n"
"    -a, --min_read_length   min read length. (not required)\n"
"    -b, --max_read_length   max read length. (not required)\n"
"\n"
"For example,\n"
"longreadqc filterfq -i input.fastq -p ./output_prefix \n"
"longreadqc filterfq -l input.fastqs.list -p ./output_prefix -n 10 \n"
);

    return 0;

}

static inline int is_bad_fastq (char **lines, int32_t min_read_len, int32_t max_read_len)
{
    int seq_len; 
    int qual_len;
    if (lines[0][0] != '@') { return 1; }
    if (lines[2][0] != '+') { return 2; }

    seq_len  = strlen(lines[1])-1;  // the last char is '\n'
    qual_len = strlen(lines[3])-1;  // the last char is '\n'

    if (seq_len != qual_len) { return 3; }

    if (seq_len < min_read_len ) { return 4; }
    if (seq_len > max_read_len ) { return 5; }

    return 0;
}

int filter1fastq(char * input_file, int num_split_file, FILE ** out_file_fps, int32_t min_read_len, int32_t max_read_len, std::unordered_map <std::string, int> read_id_map)
{
    char * read_seq;
    char ** lines;
    int32_t read_len;
    int ret;
    int read_idx;
    int out_file_idx;
    int error_code;
    int k;
    int num_skipped_reads;
    int num_good_reads;
    int num_dup_reads;

    static int total_num_skipped_reads = 0;
    static int total_num_good_reads = 0;
    static int total_num_dup_reads = 0;

	gzFile input_fp;
    std::string read_id;
    
    lines = (char **) calloc (4, sizeof(char *));
    for (int i = 0; i < 4; i++)
    {
        lines[i] = (char *) calloc(MAX_READ_LENGTH, sizeof(char));
    }

    input_fp = gzopen(input_file, "r");
    if(! input_fp) {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s", input_file);
        exit(1);
    }

    read_idx = 0;
    k = 0;
    num_skipped_reads = 0;
    num_good_reads = 0;
    num_dup_reads = 0;

    while ( gzgets(input_fp, lines[k], MAX_READ_LENGTH) != NULL )
    {
        ret = 1;
        for (int i = k+1; i < 4; i++) {
            if (gzgets(input_fp, lines[i], MAX_READ_LENGTH) == NULL) { ret = 0; }
        }
        if (ret == 0) { break; }
        error_code = is_bad_fastq(lines, min_read_len, max_read_len);
        if (error_code == 0){ // good fastq
            read_id = lines[0];
            std::size_t end_pos = read_id.find_first_of(" \t\n\r");
            if (end_pos == std::string::npos){
                read_id = read_id.substr(1);
            }else{
                read_id = read_id.substr(1, end_pos-1);
            }
            if (read_id_map.count(read_id) > 0){
                num_dup_reads += 1;
                fprintf(stderr, "WARNING! Skipped duplicated reads: %s\n", read_id.c_str()); 
                continue;
            }else{
                read_id_map[read_id] = 1;
            }
            out_file_idx = read_idx % num_split_file;  
            for (int j = 0; j < 4; j++){
                fprintf(out_file_fps[out_file_idx], "%s", lines[j]);
            }
            k = 0;
            read_idx += 1;
            num_good_reads += 1;
        } else {
            num_skipped_reads += 1;
            int at_line_num = 0; 
            for (int j = 3; j > 0; j--)
            {
                if (lines[j][0] == '@'){
                    at_line_num = j;
                    break;
                }
            }
            if (at_line_num == 0) {
                k = 0;
            }else{
                for (int i = at_line_num; i < 4; i++)
                {
                    strcpy(lines[i-at_line_num], lines[i]);
                }
                k = 4-at_line_num;  
            }
        }
    }


    total_num_skipped_reads += num_skipped_reads;
    total_num_good_reads += num_good_reads;
    total_num_dup_reads += num_dup_reads;

    fprintf(stderr, "number of clean reads of this fastq: %d\n", num_good_reads);
    fprintf(stderr, "number of filtered reads of this fastq: %d\n", num_skipped_reads);
    fprintf(stderr, "number of duplicated reads of this fastq: %d\n", num_dup_reads);
    fprintf(stderr, "total number of clean reads: %d\n", total_num_good_reads);
    fprintf(stderr, "total number of filtered reads: %d\n", total_num_skipped_reads);
    fprintf(stderr, "total number of duplicated reads: %d\n", total_num_dup_reads);
    fprintf(stderr, "\n");

    gzclose(input_fp);
    
    for (int i = 0; i < 4; i++) {
        free(lines[i]);
    }

    free(lines); 

    return 0;
}

static int filter_fastqs (STRING_LIST * input_file_list, const char * out_prefix, int num_split_file, int32_t min_read_len, int32_t max_read_len)
{
    char * input_file = NULL;
    char ** out_files;
    FILE ** out_file_fps;
    std::unordered_map <std::string, int> read_id_map;

    for (int i = 0; i < 256; i++)
    {
        if (i >= 'a' && i <= 'z' || i >= 'A' && i <= 'Z' || i == '\n') 
        {
            is_letter[i] = 1;
        }
        if (i >= 33 && i <= 126 )
        {
            is_qual_value[i] = 1;
        }
    }

    out_files = (char **) calloc (num_split_file, sizeof(char *));
    if (num_split_file > 1){
        for (int i = 0; i < num_split_file; i++) {
            out_files[i] = (char *) calloc(MAX_PATH_LENGTH, sizeof(char)); 
            sprintf(out_files[i], "%s.clean_%d.fastq", out_prefix, i);
        }
    }else if (num_split_file == 1){
        out_files[0] = (char *) calloc(MAX_PATH_LENGTH, sizeof(char)); 
        sprintf(out_files[0], "%s.clean.fastq", out_prefix);
    }

    out_file_fps = (FILE **) calloc(num_split_file, sizeof(FILE *));
    for (int i = 0; i < num_split_file; i++) {
        out_file_fps[i] = fopen(out_files[i], "w");
        if (out_file_fps[i] == NULL){
            fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", out_files[i]);
            exit(1);
        }
    }

    for (int i = 0; i < input_file_list->size; i++) {
        input_file = input_file_list->data_list[i];
        fprintf(stderr, "processing input file: %s\n", input_file);
        filter1fastq(input_file, num_split_file, out_file_fps, min_read_len, max_read_len, read_id_map); 
    }

    // free memory  
    for (int i = 0; i < num_split_file; i++) {
        free( out_files[i] );
        fclose(out_file_fps[i]);
    }

    free(out_files);
    free(out_file_fps);

    return 0;
}

int main_filterfq (int argc, char * argv[])
{
    char * input_file = NULL; 
    char * input_list_file  = NULL;
    char * out_prefix = NULL;
    char * line = NULL;
    int32_t min_read_len = 0;
    int32_t max_read_len = INT32_MAX-2;
    int     num_split_file = 1;

    STRING_LIST * input_file_list;

    line = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));

    if (argc < 5) { 
        filterfq_usage(); 
        return 1; 
    }

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 ) {
        filterfq_usage();
        return 0;
    }

    input_file_list = init_string_list(1000, MAX_PATH_LENGTH);

    int index = 1;
    while (index < argc){

        if (strcmp(argv[index], "--input_file") == 0 || strcmp(argv[index], "-i" ) == 0){
            if (index +1 >= argc){ filterfq_usage(); return 1;};
            input_file = argv[index+1];
            index += 2;
        }else if (strcmp(argv[index], "--input_list_file") == 0 || strcmp(argv[index], "-l" ) == 0){
            if (index +1 >= argc){ filterfq_usage(); return 1;};
            input_list_file = argv[index+1];
            index += 2;
        }else if (strcmp(argv[index], "--out_prefix") == 0 || strcmp(argv[index], "-p" ) == 0){
            if (index +1 >= argc){ filterfq_usage(); return 1;};
            out_prefix = argv[index+1];
            index += 2;
        }else if (strcmp(argv[index], "--num_split_file") == 0 || strcmp(argv[index], "-n" ) == 0){
            if (index +1 >= argc){ filterfq_usage(); return 1;};
            num_split_file = atoi(argv[index+1]);
            index += 2;
        }else if (strcmp(argv[index], "--min_read_length") == 0 || strcmp(argv[index], "-a" ) == 0){
            if (index +1 >= argc){ filterfq_usage(); return 1;};
            min_read_len = atoi(argv[index+1]);
            index += 2;
        }else if (strcmp(argv[index], "--max_read_length") == 0 || strcmp(argv[index], "-b" ) == 0){
            if (index +1 >= argc){ filterfq_usage(); return 1;};
            max_read_len = atoi(argv[index+1]);
            index += 2;
        }else{
            filterfq_usage(); 
            return 1;
        }
    }

    if (input_file != NULL && input_list_file != NULL){
        fprintf(stderr, "--input_file and --input_list_file cannot be specified together");
        filterfq_usage();
        exit(1);
    }

    if (input_file == NULL && input_list_file == NULL){
        fprintf(stderr, "ERROR! No input files. Both --input_file and --input_list_file were not specified.");
        filterfq_usage();
        exit(1);
    }

    if (out_prefix == NULL) { 
        fprintf(stderr, "ERROR! --out_prefix was not specified.");
        filterfq_usage();
        exit(1);
    }
    if (num_split_file < 1){
        fprintf(stderr, "ERROR! --num_split_file should be a positive integer.");
        filterfq_usage();
        exit(1);
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
    /*
    for (i = 0; i < input_file_list->size; i++)
    {
        fprintf(stdout, "%s\n", input_file_list->data_list[i]);
    }
    */
    fprintf(stdout, "output prefix: %s\n", out_prefix);

    filter_fastqs (input_file_list, out_prefix, num_split_file, min_read_len, max_read_len);

    return 0;

}
