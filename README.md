# NextSV3: an automated pipeline for structrual variation detection from long-read sequencing. 

NextSV3 uses Minimap2 to do read mapping and uses two state-of-the-art SV callers (Sniffles and cuteSV) to do SV calling.

## Installation

### Prerequisites

- Python 3.8
- [samtools](https://github.com/samtools/samtools) (v1.9 or later)
- [minimap2](https://github.com/lh3/minimap2) (v2.22 or later)
- [sniffles](https://github.com/fritzsedlazeck/Sniffles) (v2.0.7 or later)
- [cuteSV](https://github.com/tjiangHIT/cuteSV) (v2.0.1 or later)
- [GCC](https://gcc.gnu.org/) (v4.8.2 or later)

Conda environment is not required but we recommend you install the prerequisites in a new conda environment to avoid potential dependency issues. 

```
conda create -n nextsv3 python=3.8
conda activate nextsv3
conda install -c bioconda sniffles=2.0.7 cutesv=2.0.1 samtools=1.16 minimap2=2.24
```

### Installation of NextSV3
```
git clone --recursive   https://github.com/Nextomics/nextsv.git
cd nextsv/
sh build.sh
```

## Usage

Qukck Example: 
```
python path/to/nextsv3.py  -i path/to/fastq/ -o ./nextsv_output -s sample_name -r path/to/hg38.fasta -t 8
```

Full Usage:
```
usage: nextsv3.py [-h] -i path/to/input_dir -o path/to/output_dir -s sample_name -r ref.fasta [-t INT] [--samtools path/to/samtools] [--minimap2 path/to/minimap2] [--sniffles path/to/sniffles] [--cutesv path/to/cuteSV]
                  [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i path/to/input_dir, --in_dir path/to/input_dir
                        (required) path to input FASTQ/FASTA directory (all fastq/fasta files in this directory will be used as input files)
  -o path/to/output_dir, --out_dir path/to/output_dir
                        (required) path to output fastq directory
  -s sample_name, --sample_name sample_name
                        (required) a unique name or id for the input sample
  -r ref.fasta, --ref_fasta ref.fasta
                        (required) path to reference genome sequence in FASTA format
  -t INT, --threads INT
                        (optional) number of threads (default: 4)
  --samtools path/to/samtools
                        (optional) path to samtools (default: using environment default)
  --minimap2 path/to/minimap2
                        (optional) path to minimap2 (default: using environment default)
  --sniffles path/to/sniffles
                        (optional) path to sniffles (default: using environment default)
  --cutesv path/to/cuteSV
                        (optional) path to cuteSV (default: using environment default)
  -v, --version         show program's version number and exit

```

## Output files

Please go to the output directory. There will be a `work.sh`. Please run the `work.sh` locally or submit it to the cluster. 

## Contact

For questions/bugs/issues, please post on [GitHub](https://github.com/Nextomics/nextsv). In general, please do NOT send questions to our email. Your question may be very likely to help other users.

## Citation

Fang L, Hu J, Wang D, Wang K. [NextSV: a meta-caller for structural variants from low-coverage long-read sequencing data](https://doi.org/10.1186/s12859-018-2207-1).  Bioinformatics (2018) 19:180. DOI: 10.1186/s12859-018-2207-1


## More information

* [NextSV Homepage](https://github.com/Nextomics/nextsv)

## Copyright

NextSV is freely available for academic use. It is provided without warranty of any kind, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. No liability for the software usage is assumed.

For commercial use please contact Grandomics Biosciences (support@grandomics.com) for licensing options. 

Redistribution is allowed. Modification is allowed. Redistribution of modified version is not allowed, but users can submit a push request to github, and after reviewing the modification, we may accept it in the master branch.
