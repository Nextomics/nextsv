# include <stdio.h>
# include <stdlib.h>
# include <zlib.h>
# include <string.h>
# include <ctype.h>
# include "kseq.h"
# include "tk.h"



char * str_upper(const char * in_string) 
// return a new string in which all the letters are in upper case
{
    int l;
    int i;
    l = strlen(in_string);

    char * new_str;

    new_str = (char *) calloc(l+1, sizeof(char));

    for (i = 0; i < l; i++)
    {
        new_str[i] = toupper(in_string[i]);
    }
    return new_str;
}


STRING_LIST * init_string_list(int capacity, size_t max_string_length)
{
    int i;

    STRING_LIST * new_string_list;
    new_string_list = (STRING_LIST *) calloc(1, sizeof(STRING_LIST));
    new_string_list->size = 0;
    new_string_list->capacity = capacity;
    new_string_list->max_string_length = max_string_length;
    new_string_list->data_list = (char **)calloc(new_string_list->capacity, sizeof(char *)); 

    if (new_string_list->data_list == NULL){
        fprintf(stderr, "ERROR! Failed to alloc memory for string list (capacity=%d, max_string_length=%d)\n", capacity, max_string_length);
        exit(1);
    }
    for (i = 0; i < new_string_list->capacity; i++){
        new_string_list->data_list[i] = (char *)calloc(max_string_length, sizeof(char));
        if (new_string_list->data_list[i] == NULL){
            fprintf(stderr, "ERROR! Failed to alloc memory for string list (capacity=%d, max_string_length=%d)\n", capacity, max_string_length);
            exit(1);
        }
    }

    return new_string_list;

}

int append_string_list (STRING_LIST * string_list, char * new_string)
{
    int l;
    int i;
    l = strlen(new_string);

    if (l > string_list->max_string_length-1){
        fprintf(stderr, "ERROR! Cannot append new string to string list! Too long\n");
        exit(1);
    }
    if (string_list->capacity == string_list->size)
    {
        string_list->capacity  = string_list->capacity * 2;
        string_list->data_list = (char **) realloc (string_list->data_list, string_list->capacity*sizeof(char*));
        if (string_list->data_list[i] == NULL){
            fprintf(stderr, "ERROR! Failed to realloc memory for string list (capacity=%d, max_string_length=%d)\n", string_list->capacity, string_list->max_string_length);
            exit(1);
        }
        for (i = string_list->size; i < string_list->capacity; i++)
        {
            string_list->data_list[i] = (char *)calloc(string_list->max_string_length, sizeof(char));
            if (string_list->data_list[i] == NULL){
                fprintf(stderr, "ERROR! Failed to realloc memory for string list (capacity=%d, max_string_length=%d)\n", string_list->capacity, string_list->max_string_length);
                exit(1);
            }
        }
    }
    strcpy(string_list->data_list[string_list->size], new_string);
    string_list->size += 1;

    return 0;

}

int add1input_file(STRING_LIST * input_file_list, char * input_file)
{

    int l;
    int i;
    l = strlen(input_file);
    for (i = 0; i < l; i++)
    {   
        if (input_file[i] == '\n'){
            input_file[i] = '\0';
        }   
    }   
    append_string_list(input_file_list, input_file);

    return 0;
}

int get_read_length_statistics (int * read_length_count, int read_length_limit, int * p_n50_read_length, int * p_mean_read_length, int *p_median_read_length, int *p_n05_read_length, int * p_n95_read_length, int *p_longest_read_length)
{
    double sum_read_length;
    int i;
    double curr_sum;
    double total_num_reads;

    sum_read_length = 0;
    total_num_reads = 0;

    for (i = read_length_limit; i >= 0; i--) {
        sum_read_length += read_length_count[i] * i;
        total_num_reads += read_length_count[i];
    }

    curr_sum = 0;

    * p_n05_read_length = -1;
    * p_n50_read_length = -1;
    * p_n95_read_length = -1;
    * p_longest_read_length = -1;
    for (i = read_length_limit; i >= 0; i--)
    {   
        if (*p_longest_read_length < 0 && read_length_count[i] > 0){
            *p_longest_read_length = i;
        }
        curr_sum += read_length_count[i] * i;  
        if (*p_n05_read_length < 0 && curr_sum  >= 0.05 * sum_read_length){
            * p_n05_read_length = i;
        }
        if (*p_n50_read_length < 0 && curr_sum  >= 0.5 * sum_read_length) {   
            *p_n50_read_length = i; 
        } 
        if (*p_n95_read_length < 0 && curr_sum  >= 0.95 * sum_read_length){
            *p_n95_read_length = i;
            break;
        }
    }   

    *p_mean_read_length = sum_read_length / total_num_reads; 

    int curr_read_cnt_sum = 0;
    for (i = read_length_limit; i >= 0; i--)
    {
        curr_read_cnt_sum += read_length_count[i];
        if (curr_read_cnt_sum >= 0.5 * total_num_reads)
        {
            *p_median_read_length = i;
            break;
        }
    }

    return 0;
}

int output_read_length_histo(int * read_length_count, int read_length_limit, char * read_length_histo_file, int bin_size )
{
    int i;
    FILE * read_length_histo_fp;
    int n_bin;
    int64_t * bin_read_count;
    int64_t * bin_base_count;
    int bin_idx; 

    read_length_histo_fp = fopen(read_length_histo_file, "w");
    fprintf(read_length_histo_fp, "#read_length\tnum_reads\n");

    n_bin = read_length_limit/bin_size + 2;
    bin_read_count = (int64_t *)calloc(n_bin, sizeof(int64_t));
    bin_base_count = (int64_t *)calloc(n_bin, sizeof(int64_t));

    for (bin_idx = 0; bin_idx < n_bin; bin_idx++) { 
        bin_read_count[bin_idx] = 0;
        bin_base_count[bin_idx] = 0;
    }

    for (i = 0; i<= read_length_limit; i++)
    {   
        bin_idx = (i + bin_size/2)/bin_size;
        bin_read_count[bin_idx] += read_length_count[i]; 
        bin_base_count[bin_idx] += read_length_count[i] * i; 
    }   
    for (bin_idx = 0; bin_idx < n_bin; bin_idx++) {
        fprintf(read_length_histo_fp, "%d\t%lld\t%lld\n", bin_idx * bin_size, bin_read_count[bin_idx], bin_base_count[bin_idx]);
    }

    free(bin_read_count);
    free(bin_base_count);
    fclose(read_length_histo_fp);

    return 0;
}


char * join_path(const char *path1, const char *path2)
{
    char * joined_path = (char *) calloc(MAX_PATH_LENGTH, sizeof(char));
    sprintf(joined_path, "%s/%s", path1, path2);
    return joined_path;
}
