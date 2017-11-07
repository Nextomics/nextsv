# NextSV: a meta-caller for structrual variation from low coverage long-read sequencing.

NextSV, a meta SV caller and a computational pipeline to perform SV calling from low coverage long-read sequencing data. NextSV integrates three aligners and three SV callers and generates two integrated call sets (sensitive/stringent) for different analysis purpose. The output of NextSV is in ANNOVAR-compatible bed format. Users can easily perform downstream annotation using ANNOVAR and disease gene discovery using Phenolyzer.

## Features

* Suitable for identifying structural variants from Mendelian diseases.

* Fast and easy customization

## Supported aligners

BLASR, BWA-MEM, NGMLR

## Supported SV callers

PBHoney, Sniffles

## Installation

Prerequisites:
   
zlib-dev, cmake, gcc/g++(>=4.8.2), pip, bwa, samtools

Downloading NextSV: 

```
git clone https://github.com/Nextomics/nextsv.git
```

Installation of PBHoney:
```
cd nextsv/
pip install --user argparse
pip install --user pysam
pip install --user networkx
pip install --user h5py
pip install --user PyIntervalTree
sh install_pbhoney.sh
```
Installation of Sniffles:
```
sh install_sniffles.sh
```
If you have sniffles installed on your system, you can skip this step and specify the path to sniffles in the config file.

For user's convenience, we provided compiled binary files of NGMLR and BLASR. The binary file of NGMLR was released by its author (https://github.com/philres/ngmlr). The binary file of BLASR was compiled by us. The copyright of NGMLR and BLASR can be found in the Copyright section of this file. We tested the two binary files on CentOS 6.5 and Ubuntu 14.


## Usage
```
python nextsv.py [config_file]
```
A template config file can be found in example.config. The following parameters can be set in the config file:

`sample_name`: sample name. It will be a part of the prefix of the output files. 

`input_file_list`: path to a file that contains names of input files. Each line contains one input file. Examples of input file list can be found in `example.fastq.fofn`. 

`out_dir`: full path to output directory

`n_thread`: number of threads (CPU cores)

`mode`: running mode of NextSV. NextSV currently supports two modes: multiprocessing or sge. `mode=multiprocessing` means that NextSV will run in parallel using multiple cores (equal to `n_thread`) on a single machine. `mode=sge` means that NextSV will submit the jobs to the SGE (Sun Grid Engine) cluster. If `mode=sge` is specified, users need to provide the submission command with parameters (`job_submission_command`).

`job_submission_command`: job submission command with parameters. Please specify number of threads in the parameters. For example, `job_submission_command=[qsub -V -cwd -S /bin/bash -l h_vmem=4G -pe smp 12]` means that the job will be submitted using qsub and the parameters for qsub is "-V -cwd -S /bin/bash -l h_vmem=4G -pe smp 12", which means the job will use 12 threads and 4G memory for each thread. NextSV will copy this command while submitting the shell script file. 

`enable_PBHoney_Spots`: whether to enable SV calling using PBHoney-Spots (1 for enable, 0 for disable)

`spots_threshold`: threshold parameter of PBHoney-Spots (default=2)

`spots_minErrReads`: minimal supporting reads for PBHoney-Spots (default=2)

`spots_consensus`: method for polishing consensus for PBHoney-Spots(pbdagcon or None, default=None)


`enable_PBHoney_Tails`: whether to enable SV calling using PBHoney-Tails (1 for enable, 0 for disable)

`tails_minBreads`: minimum supporting reads for PBHoney-Tails (default=2)

`tails_minZMWs`: minimum number of unique ZMWs for PBHoney-Tails (default=2)

`tails_buffer`: buffer parameter of PBHoney-Tails(default=600)


`enable_bwa_Sniffles`: whether to enable SV calling using the combination of BWA and Sniffles (1 for enable, 0 for disable)

`enable_ngmlr_Sniffles`: whether to enable SV calling using the combination of NGMLR and Sniffles (1 for enable, 0 for disable)

`sniffles_min_support`: minimum supporting reads for Sniffles (default=2)

`sniffles_max_distance`:maximum distance to group SV together by Sniffles (default=600)

`bwa`: full path to bwa binary file, required

`ngmlr`: full path to ngmlr binary file

`sniffles`: full path to sniffles binary file

`bash5tools`: full path to installed bash5tools.py on your system, required if your input files are in hdf5 format. This is for extracting fastq from hdf5 files. NextSV currently support pbh5tools (bash5tools.py) version 0.8.0

`minLength`: minimal subread length for extracting fastq from hdf5 file (default: 500). 

`minReadScore`: minReadScore for extracting fastq from hdf5 file (default: 0.75).

`samtools`: full path to samtools binary file, , required

`ref_fasta`: full path to reference genome fasta file (for BWA and NGMLR aligners. If `enable_bwa_Sniffles=1` is specified, the fasta file should be pre-indexed by BWA. If `enable_bwa_Sniffles=0` is specified, the fasta file does not need to be pre-indexed)

`ref_blasr`: full path to reference genome fasta file (for blasr aligner)

`ref_sa_blasr`: full path to precomputed suffix array of reference genome fasta file for blasr aligner (for information of generating the suffix array, please refer to https://github.com/PacificBiosciences/blasr/blob/master/README.MANUAL.md)

Please use full paths in the config file and input list file.

## Contact

For questions/bugs/issues, please post on [GitHub](https://github.com/Nextomics/nextsv). In general, please do NOT send questions to our email. Your question may be very likely to help other users.

## Citation

Fang L, Hu J, Wang D, Wang K. [NextSV: a computational pipeline for structural variation analysis from low-coverage long-read sequencing](http://www.biorxiv.org/content/early/2017/07/17/092544). bioRxiv, doi: 10.1101/092544

## More information

* [NextSV Homepage](https://github.com/Nextomics/nextsv)

## Copyright

### Copyright of NextSV
NextSV is freely available for academic use. It is provided without warranty of any kind, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. No liability for the software usage is assumed.

For commercial use please contact Grandomics Biosciences (support@grandomics.com) for licensing options. 

Redistribution is allowed. Modification is allowed. Redistribution of modified version is not allowed, but users can submit a push request to github, and after reviewing the modification, we may accept it in the master branch.

### Copyright of NGMLR
The MIT License (MIT)

Copyright (c) 2017 Philipp Rescheneder

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### Copyright of BLASR
Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.

All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted (subject to the limitations in the
disclaimer below) provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

 * Neither the name of Pacific Biosciences nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.


