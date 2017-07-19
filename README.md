# NextSV: a computational pipeline for structrual variation detection from low coverage long-read sequencing.

NextSV, a meta SV caller and a computational pipeline to perform SV calling from low coverage long-read sequencing data. NextSV integrates three aligners and three SV callers and generates two integrated call sets (sensitive/stringent) for different analysis purpose. The output of NextSV is in ANNOVAR-compatible bed format. Users can easily perform downstream annotation using ANNOVAR and disease gene discovery using Phenolyzer.

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
   
   zlib-dev, cmake, gcc/g++(>=4.8.2), pip

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

