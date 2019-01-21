# NextSV2: an automated pipeline for structrual variation detection from long-read sequencing.

NextSV2 uses two aligners (Minimap2 and NGMLR) to do read mapping and uses Sniffles to do SV calling.

## Installation

Prerequisites:
sniffles, samtools (version 1.3 or later), gcc/g++(>=4.8.2), python 2.7

Installation:

```
git clone --recursive   https://github.com/Nextomics/nextsv.git
cd nextsv/
sh install.sh
```

## Usage

```
python nextsv2.py <config> <sample_name> <input_fastq_folder>
```

`sample_name` is a unique name or id for the input sample. 

`input_fastq_folder` is the folder that contains all input fastq files. NextSV2 will search the folder and all fastq files (.fastq or .fq) found in this folder will be treated as input files. 

A template config file can be found in example.config. 

Please use full paths in the config file and input list file.

## Output files

Please go to the output directory. There will be a `work.sh`. Please run the `work.sh` locally or submit it to the cluster. 
The output VCF files is in the `sniffles_calls` directory.

## Contact

For questions/bugs/issues, please post on [GitHub](https://github.com/Nextomics/nextsv). In general, please do NOT send questions to our email. Your question may be very likely to help other users.

## Citation

Fang L, Hu J, Wang D, Wang K. [NextSV: a meta-caller for structural variants from low-coverage long-read sequencing data](https://doi.org/10.1186/s12859-018-2207-1).  Bioinformatics (2018) 19:180. DOI: 10.1186/s12859-018-2207-1


## More information

* [NextSV Homepage](https://github.com/Nextomics/nextsv)

## Copyright

### Copyright of NextSV
NextSV is freely available for academic use. It is provided without warranty of any kind, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. No liability for the software usage is assumed.

For commercial use please contact Grandomics Biosciences (support@grandomics.com) for licensing options. 

Redistribution is allowed. Modification is allowed. Redistribution of modified version is not allowed, but users can submit a push request to github, and after reviewing the modification, we may accept it in the master branch.


