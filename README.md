# NextSV3: automated structrual variation detection from long-read sequencing using state-of-the-art tools. 

NextSV3 uses Minimap2 to do read mapping and uses two state-of-the-art SV callers (Sniffles and cuteSV) to do SV calling.

## Installation

### Prerequisites

- Python 3.8
- [samtools](https://github.com/samtools/samtools) (v1.9 or later)
- [minimap2](https://github.com/lh3/minimap2) (v2.22 or later)
- [sniffles](https://github.com/fritzsedlazeck/Sniffles) (v2.0.7 or later)
- [cuteSV](https://github.com/tjiangHIT/cuteSV) (v2.0.1 or later)
- [GCC](https://gcc.gnu.org/) (v4.8.2 or later)

Conda environment is not required but we highly recommend you install the prerequisites in a new conda environment to avoid potential dependency issues. 

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

### Quick Example

```
# Oxford Nanopore reads
python path/to/nextsv3.py  -i path/to/fastq/folder/ -o ./nextsv_output -s unique_sample_name -r path/to/hg38.fasta -t 8 -p ont -e nextsv3

# PacBio CLR reads
python path/to/nextsv3.py  -i path/to/fastq/folder/ -o ./nextsv_output -s unique_sample_name -r path/to/hg38.fasta -t 8 -p clr -e nextsv3

# PacBio HiFi reads
python path/to/nextsv3.py  -i path/to/fastq/folder/ -o ./nextsv_output -s unique_sample_name -r path/to/hg38.fasta -t 8 -p hifi -e nextsv3
```

Memory consumption of the pipeline depends on number of threads and the size of the reference genome. For human genomes (3Gb), we recommend 4GB memory per thread. 

### Full Usage
```
usage: nextsv3.py [-h] -i path/to/input_dir -o path/to/output_dir -s sample_name -r ref.fasta -p sequencing_platform [-t INT] [-e conda_env] [--samtools path/to/samtools] [--minimap2 path/to/minimap2]
                  [--sniffles path/to/sniffles] [--cuteSV path/to/cuteSV] [-v]

nextsv3: an automated pipeline for structrual variation detection from long-read sequencing. Contact: Li Fang(fangli2718@gmail.com)

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
  -p sequencing_platform, --platform sequencing_platform
                        (required) sequencing platform. Three valid values: ont/clr/hifi ont: oxford nanopore; clr: PacBio CLR; hifi: PacBio HiFi/CCS
  -t INT, --threads INT
                        (optional) number of threads (default: 4)
  -e conda_env, --conda_env conda_env
                        (optional) conda environment name (default: NULL)
  --samtools path/to/samtools
                        (optional) path to samtools (default: using environment default)
  --minimap2 path/to/minimap2
                        (optional) path to minimap2 (default: using environment default)
  --sniffles path/to/sniffles
                        (optional) path to sniffles (default: using environment default)
  --cuteSV path/to/cuteSV
                        (optional) path to cuteSV (default: using environment default)
  -v, --version         show program's version number and exit
```

## Output files

SV calls of sniffles and cuteSV will be generated in the `out_dir/3_SV_calls` folder.

## Contact

If you have any questions/suggestions, please feel free to post it on the [Issue page](https://github.com/Nextomics/nextsv/issues). You may also email me at `fangli2718@gmail.com`. 

## Citation

Fang L, Hu J, Wang D, Wang K. [NextSV: a meta-caller for structural variants from low-coverage long-read sequencing data](https://doi.org/10.1186/s12859-018-2207-1).  Bioinformatics (2018) 19:180. DOI: 10.1186/s12859-018-2207-1

