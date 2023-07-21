# LongReadQC
LongReadQC is a quality control tool for long reads. LongReadQC is written in C. 
LongReadQC is in development

## Usage
LongReadQC currently has two functions: 

1) QC for fastq file. LongReadQC will generate a text file showing basic information of the input reads. 

```         
Usage:   longreadqc fq [options]  
Options:
    -h, --help              output this usage information
    -i, --input_file        path to the input file, which can be fasta, fastq or bam. Use this argument if only you have only one input file
    -l, --input_list_file   a file that contains the paths of all the input files. Usage this argument if you have multiple input files
    -p, --out_prefix     prefix of output files. (required)
    -d, --out_dir        output directory. (default: ./longreadqc_out/)

For example,
longreadqc fq -i FileName.fastq -d ./FileName_longreadqc_out/ -p sample01 
longreadqc fq -l SampleName.fastqs.list -d ./SampleName_longreadqc_out/ -p sample02 
```  

2) Filtering reads. LongReadQC will remove reads that are illegal, for example, reads that have different length in base and quality. Users can also set min_read_length and max_read_length to get the reads in the intended length. 

```
Usage:   longreadqc filterfq [options]  
Options:
    -h, --help              output this usage information
    -i, --input_file        path to the input fastq file. Use this argument if only you have only one input file
    -l, --input_list_file   a file that contains the paths of all the input fastq files. Usage this argument if you have multiple input files
    -p, --out_prefix        prefix of the output file. (required)
    -n, --num_split_file    split the output file to n files. (not required, default: n=1)
    -a, --min_read_length   min read length. (not required)
    -b, --max_read_length   max read length. (not required)

For example,
longreadqc filterfq -i input.fastq -p ./output_prefix 
longreadqc filterfq -l input.fastqs.list -p ./output_prefix -n 10 
```

## Installation

```
git clone https://github.com/fangli08/longreadqc.git
cd longreadqc/
make
```

You can copy the `longreadqc` binary file to the folder in your PATH.


