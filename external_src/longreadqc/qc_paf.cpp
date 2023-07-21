# include <stdio.h>
# include <stdlib.h>
# include <zlib.h>
# include <string.h>
# include <ctype.h>
# include <stdint.h>

# include <sys/types.h>
# include <sys/stat.h>

# include "kseq.h"
# include "tk.h"

# define MAX_INDEL_LENGTH 10000

static int paf_usage()
{
    printf(
"\n"
"Program: longreadqc (Quality control tool for long read sequencing)\n\n"
"Usage:   longreadqc paf \n"
"Options:\n"
"    -h, --help              output this usage information\n"
"    -i, --input_file        path to the input file, which should be a paf file\n"
"    -p, --out_prefix        prefix of the output files. (default: ./InputFileName.statistics.txt)\n");

    return 0;

}

int get_var_stat_from_short_cs_string(const char * cs_string, int64_t * ins_length_list, int64_t * del_length_list, int64_t max_indel_length, int64_t * p_n_match, int64_t * p_n_mismatch, int64_t *p_n_ins_event, int64_t * p_n_del_event)
{

    int opera_table[256] = {0};
    int64_t n_match;
    int64_t l_ins, l_del;
    int64_t i, j;
    int64_t cs_len;
    char n_match_string[32] = {0};

    opera_table[':'] = 1;
    opera_table['+'] = 1;
    opera_table['-'] = 1;
    opera_table['*'] = 1;
    opera_table['~'] = 1;

    i = 0;
    cs_len = strlen(cs_string);
    while (i < cs_len)
    {
        if (cs_string[i] == '*'){
            j = i + 1;
            while (j < cs_len){
                if (opera_table[cs_string[j]] == 1) { break;}
                j++;
            } // eventually j equals to cs_len or a position of :+-* 

            (*p_n_mismatch)++;
            i = j;
        }else if (cs_string[i] == '+'){
            (*p_n_ins_event)++;

            j = i + 1;
            while (j < cs_len){
                if (opera_table[cs_string[j]] == 1) { break;}
                j++;
            } // eventually j equals to cs_len or a position of :+-* 

            l_ins = j - i - 1;
            if (l_ins <= max_indel_length-1) {
                ins_length_list[l_ins]++; 
            }else{
                ins_length_list[max_indel_length-1]++; 
            }
            i = j;
        }else if (cs_string[i] == '-'){
            (*p_n_del_event)++;

            j = i + 1;
            while (j < cs_len){
                if (opera_table[cs_string[j]] == 1) { break;}
                j++;
            } // eventually j equals to cs_len or a position of :+-* 

            l_del = j - i - 1;
            if (l_del <= max_indel_length-1) {
                del_length_list[l_del]++; 
            }else{
                del_length_list[max_indel_length-1]++; 
            }
            i = j;
        }else if (cs_string[i] == ':'){
            j = i + 1;
            while (j < cs_len){
                if (opera_table[cs_string[j]] == 1) { break;}
                j++;
            } // eventually j equals to cs_len or a position of :+-* 

            memcpy(n_match_string, cs_string + i + 1, j - i - 1);
            n_match_string[j-i-1] = 0;
            n_match = atoi(n_match_string);
            *p_n_match += n_match;
            //fprintf(stderr, "i=%d,j=%d,n_match_string=%s,n_match=%d, total_n_match=%d\n", i, j, n_match_string, n_match, *p_n_match);
            i = j;
        }else{
            i++;
        }
    }

    return 0;

}

int qc_paf(const char * input_file, const char * out_prefix)
{
    gzFile input_fp;

    char * line;
    char * q_name, * t_name;
    int64_t q_len, t_len, q_start, q_end, t_start, t_end;
    char strand[10];
    int64_t n_residue_matches, alignment_block_length, map_qual;
    int64_t n_match, n_mismatch, n_ins_event, n_del_event;
    int64_t line_len = 0;
    int64_t * ins_length_list, * del_length_list;
    int64_t min_mapq = 30;
    int64_t start_idx, end_idx;
    char * out_file;
    char * cs_string;
    FILE * out_fp;
    int64_t total_ins_len, total_del_len;
    int64_t total_50_ins_len, total_50_del_len;
    int n_skipped_alignments;

    ins_length_list = (int64_t *) calloc (MAX_INDEL_LENGTH, sizeof(int64_t));
    del_length_list = (int64_t *) calloc (MAX_INDEL_LENGTH, sizeof(int64_t));
    line = (char *) calloc(MAX_PAF_LENGTH, sizeof(char));
    cs_string = (char *) calloc(MAX_PAF_LENGTH, sizeof(char));
    q_name = (char *) calloc(MAX_PAF_LENGTH, sizeof(char));
    t_name = (char *) calloc(MAX_PAF_LENGTH, sizeof(char));
    out_file = (char *)calloc(MAX_PATH_LENGTH, sizeof(char));

    input_fp = gzopen(input_file, "r");
    n_match = n_mismatch = n_ins_event = n_del_event = 0;
    if(!input_fp) {
        fprintf(stderr, "ERROR! Failed to gzopen file for reading: %s", input_file);
        exit(1);
    }

    n_skipped_alignments = 0;
    while (gzgets(input_fp, line, MAX_PAF_LENGTH) != NULL )
    {
        sscanf(line, "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%*s", q_name, &q_len, &q_start, &q_end, strand, t_name, &t_len, &t_start, &t_end, &n_residue_matches, &alignment_block_length, &map_qual); 
        if (map_qual < min_mapq){ 
            n_skipped_alignments++;
            continue;
         }
        line_len = strlen(line);
        for (int i = 0; i < line_len; i++)
        {
            if (line[i] == '\t' && line[i+1] == 'c' && line[i+2] == 's' && line[i+3] == ':' && line[i+4] == 'Z' && line[i+5] == ':')
            {
                start_idx = i + 6;
                end_idx = start_idx;
                while (end_idx < line_len && line[end_idx] != '\t')
                {
                    end_idx ++; // end_idx itself is not included.
                }
                
                for (int j = start_idx; j < end_idx; j++)
                {
                    cs_string[j-start_idx] = line[j];
                }
                cs_string[end_idx-start_idx] = '\0';
                get_var_stat_from_short_cs_string(cs_string, ins_length_list, del_length_list, MAX_INDEL_LENGTH, &n_match, &n_mismatch, &n_ins_event, &n_del_event);

                break;
            }
        }
    }

    fprintf(stderr,"number of skipped alignments:%d\n", n_skipped_alignments);
    sprintf(out_file, "%s.varstat.txt", out_prefix);
    out_fp = fopen(out_file, "w");
    if (out_fp == NULL){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", out_file);
        exit(1);
    }

    total_ins_len = total_del_len = 0;
    for (int i = 1; i < MAX_INDEL_LENGTH; i++)
    {
        total_ins_len += ins_length_list[i] * i;
        total_del_len += del_length_list[i] * i;
    }

    total_50_ins_len = total_50_del_len = 0;
    for (int i = 1; i <= 50; i++)
    {
        total_50_ins_len += ins_length_list[i] * i;
        total_50_del_len += del_length_list[i] * i;
    }


    fprintf(out_fp, "number of matched bases\t%lld\n", n_match);
    fprintf(out_fp, "number of mismatched bases\t%lld\n", n_mismatch);
    fprintf(out_fp, "number of insertion events\t%lld\n", n_ins_event);
    fprintf(out_fp, "number of deletion events\t%lld\n", n_del_event);
    fprintf(out_fp, "number of insertion bases(ins len <= 50 bp)\t%lld\n", total_50_ins_len);
    fprintf(out_fp, "number of deletion bases (del len <= 50 bp)\t%lld\n", total_50_del_len);
    fprintf(out_fp, "number of insertion bases\t%lld\n", total_ins_len);
    fprintf(out_fp, "number of deletion bases\t%lld\n", total_del_len);
    fprintf(out_fp, "\n\n\n");
    fprintf(out_fp, "Length\tnum_ins_events\tnum_del_events\tnum_ins_bases\tnum_del_bases\n");
    for (int64_t i = 1; i < 50; i++) {
        fprintf(out_fp, "%lld\t%lld\t%lld\t%lld\t%lld\n", i, ins_length_list[i], del_length_list[i], ins_length_list[i]*i, del_length_list[i]*i);
    }
    
    free(line); 
    free(cs_string);
    free(ins_length_list);
    free(del_length_list);
    free(q_name);
    free(t_name);
    free(out_file);
    
    return 0;
}

int main_qcpaf(int argc, char * argv[])
{

    char * input_file = NULL;
    char * out_prefix = NULL;

    if (argc < 2) {
        paf_usage();
        return 1;
    }

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 ) { 
        paf_usage();
        return 0;
    }

    int index = 1;
    while (index < argc){

        if (strcmp(argv[index], "--input_file") == 0 || strcmp(argv[index], "-i" ) == 0){
            if (index +1 >= argc){ paf_usage(); return 1;};
            input_file = argv[index+1];
            index += 2;
        }else if (strcmp(argv[index], "--out_prefix") == 0 || strcmp(argv[index], "-p" ) == 0){
            if (index +1 >= argc){ paf_usage(); return 1;};
            out_prefix = argv[index+1];
            index += 2;
        }else{
            paf_usage();
            return 1;
        }
    }

    qc_paf(input_file, out_prefix);
    return 0;
}
