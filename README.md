# NextSV: Automated structrual variation detection for long-read sequencing.

NextSV is an computational pipeline that allows structural variant (SV) calling from PacBio sequencing data using PBhoney and Sniffles. NextSV takes FASTA or FASTQ files as input. Once the SV caller is selected by user, NextSV automatically chooses the compatible aligner and performs mapping. The alignments will be automatically sorted and then presented to the SV caller. Users can change the parameters by modifying its configuration file. When the analysis is finished, NextSV will examine the FASTA/FASTQ, BAM, and result files and generate a report showing various statistics. If more than both callers are selected, NextSV will format the raw result files (.tails, .spots, or .vcf files) into bed files and generate the intersection or union call set for the purpose of higher accuracy or sensitivity.

## Features

* Suitable for identifying structural variants from Mendelian diseases.

* Automatic installation of all aligners and SV callers.

* Fast and easy customization

* Sun Grid Engine (SGE) integration

## Supported aligners

   BLASR, BWA-MEM, NGMLR

## Supported SV callers

   PBHoney, Sniffles

## Installation

   Prerequisites:
   zlib-dev, cmake, gcc/g++(>=4.8.2)

   Download NextSV: 

   ```
   git clone https://github.com/Nextomics/nextsv.git
   ```

   Install:
   ```
   cd nextsv/
   sh install.sh
   cp ~/.bashrc  ~/.bashrc.bak
   cat setup-env.sh >> ~/.bashrc

   ```

## Usage
   ```
   perl nextsv.pl [config_file]
   ```
   for configuration file, please refer to example.config

## Contact

For questions/bugs/issues, please post on [GitHub](https://github.com/Nextomics/nextsv). In general, please do NOT send questions to our email. Your question may be very likely to help other users.

## Citation

Fang L, Hu J, Wang D, Wang K. [Evaluation on Efficient Detection of Structural Variants at Low Coverage by Long-Read Sequencing](http://biorxiv.org/content/early/2016/12/09/092544). bioRxiv, doi: 10.1101/092544

## More information

* [NextSV Homepage](https://github.com/Nextomics/nextsv)

## Copyright

NextSV is freely available for academic use. It is provided without warranty of any kind, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. No liability for the software usage is assumed.

For commercial use please contact Grandomics Biosciences (support@grandomics.com) for licensing options. 

Redistribution is allowed. Modification is allowed. Redistribution of modified version is not allowed, but users can submit a push request to github, and after reviewing the modification, we may accept it in the master branch.

