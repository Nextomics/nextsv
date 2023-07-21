#ifndef TK_H 
#define TK_H

#define MAX_READ_LENGTH 10485760
#define MAX_PATH_LENGTH 10240
#define MAX_PAF_LENGTH 10485760

typedef struct {
	char ** data_list;
	size_t size;
	size_t capacity;
	size_t max_string_length;
} STRING_LIST;

char * str_upper(const char * in_string);

STRING_LIST * init_string_list (int capacity, size_t max_string_length);
int add1input_file(STRING_LIST * input_file_list, char * input_file);
int get_read_length_statistics (int * read_length_count, int read_length_limit, int * p_n50_read_length, int * p_mean_read_length, int *p_median_read_length, int *p_n05_read_length, int * p_n95_read_length, int *p_longest_read_length);
int output_read_length_histo(int * read_length_count, int read_length_limit, char * read_length_histo_file, int bin_size);
char * join_path(const char *path1, const char *path2);

#endif 
