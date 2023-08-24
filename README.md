# NextSV3: automated structrual variation detection from long-read sequencing using state-of-the-art tools

NextSV3 uses [Minimap2](https://github.com/lh3/minimap2)/[Winnowmap](https://github.com/marbl/Winnowmap)/[NGMLR](https://github.com/philres/ngmlr) to do read mapping and uses two state-of-the-art SV callers ([Sniffles2](https://github.com/fritzsedlazeck/Sniffles) and [cuteSV2](https://github.com/tjiangHIT/cuteSV)) to do SV calling.

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Installation](#installation)
- [Usage](#usage)
  - [Aligners](#aligners)
  - [Quick Example](#quick-example)
  - [Output files](#output-files)
  - [Full Usage](#full-usage)
- [Contact](#contact)
- [Citation](#citation)
- [Included tools](#open-source-software-tools-included-in-this-repository)

## Installation

**Prerequisites**
- Operating System: Linux
- Python 3.8 or later
- [samtools](https://github.com/samtools/samtools) (v1.9 or later)
- [sniffles](https://github.com/fritzsedlazeck/Sniffles) (v2.0.7 or later)
- [cuteSV](https://github.com/tjiangHIT/cuteSV) (v2.0.1 or later)

Conda environment is not required but we highly recommend you install the prerequisites in a new conda environment to avoid potential dependency issues. 

```
conda create -n nextsv3 python=3.8
conda activate nextsv3
pip install "setuptools<58.0"
pip install sniffles
pip install cuteSV
```

**Installation of NextSV3**

If your OS is Linux, you can acquire precompiled binaries from the [release page](https://github.com/Nextomics/nextsv/releases) with:

```
wget https://github.com/Nextomics/nextsv/releases/download/v3.1.0/NextSV-3.1.0-linux-x86_64.tar.gz
tar xzf NextSV-3.1.0-linux-x86_64.tar.gz
```

If you want to compile from the source code, you need to have the following tools/libraries installed:
- GCC (v7.0 or later)
- GNU make
- CMake
- zlib development files 

```
git clone https://github.com/Nextomics/nextsv.git
cd nextsv
make
```

If you are in China, you can clone NextSV from the following mirror with much better speed:
```
git clone https://e.coding.net/fanglab/LongReads/nextsv.git
cd nextsv
make
```


## Usage

### Aligners

The default aligner is `minimap2` as it is faster. But you can use the `-a`/`--aligners` argument to specifiy aligners. 
You can use any combination of the following 3 aligners: `minimap2`, `winnowmap`, `ngmlr` (separated by the `+` sign):

```
# use minimap2 only (default) 
-a minimap2

# use minimap2 and winnowmap
-a minimap2+winnowmap

# use minimap2, winnowmap and ngmlr
-a minimap2+winnowmap+ngmlr

# use winnowmap and ngmlr
-a winnowmap+ngmlr

# use minimap2 and ngmlr
-a minimap2+ngmlr

# use ngmlr only
-a ngmlr

# use minimap2 only
-a minimap2
```

Please note that if you use more aligners, it would take more time to finish the pipeline.

### Quick Example

```
# PacBio HiFi reads
python path/to/nextsv3.py -a minimap2+winnowmap+ngmlr -i path/to/fastq/folder/ -o ./nextsv_output -s unique_sample_name -r path/to/hg38.fasta -t 8 -p hifi -e nextsv3

# PacBio CLR reads
python path/to/nextsv3.py -a minimap2+winnowmap+ngmlr -i path/to/fastq/folder/ -o ./nextsv_output -s unique_sample_name -r path/to/hg38.fasta -t 8 -p clr -e nextsv3

# Oxford Nanopore reads
python path/to/nextsv3.py -a minimap2+winnowmap+ngmlr -i path/to/fastq/folder/ -o ./nextsv_output -s unique_sample_name -r path/to/hg38.fasta -t 8 -p ont -e nextsv3
```

Memory consumption of the pipeline depends on number of threads and the size of the reference genome. For human genomes (3Gb), we recommend 4GB memory per thread. 

### Output files

NextSV will generate a `work.sh` in the output directory. Run this `work.sh` and you will get output files. SV calls of sniffles and cuteSV will be generated in the `out_dir/3_SV_calls` folder.

### Full Usage
```
usage: nextsv3.py [-h] -i path/to/input_dir -o path/to/output_dir -s sample_name -r ref.fasta -p sequencing_platform -a aligners_to_use [-t INT]
                  [-e conda_env] [--samtools path/to/samtools] [--sniffles path/to/sniffles] [--cuteSV path/to/cuteSV]
                  [--minimap2 path/to/minimap2] [-v]

nextsv3: an automated pipeline for structrual variation detection from long-read sequencing. Contact: Li Fang(fangli2718@gmail.com)

options:
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
                        (required) sequencing platform. Three valid values: ont/clr/hifi ont: oxford nanopore; clr: PacBio CLR; hifi: PacBio
                        HiFi/CCS
  -a aligners_to_use, --aligners aligners_to_use
                        (optional) which aligner(s) to use. Three supported aligners: minimap2, ngmlr, winnowmap. Use "+" to combine multiple
                        aligners. Examples: minimap2+ngmlr, winnowmap+minimap2, minimap2+ngmlr+winnowmap, minimap2, ngmlr, winnowmap
  -t INT, --threads INT
                        (optional) number of threads (default: 4)
  -e conda_env, --conda_env conda_env
                        (optional) conda environment name (default: NULL)
  --samtools path/to/samtools
                        (optional) path to samtools (default: using environment default)
  --sniffles path/to/sniffles
                        (optional) path to sniffles (default: using environment default)
  --cuteSV path/to/cuteSV
                        (optional) path to cuteSV (default: using environment default)
  --minimap2 path/to/minimap2
                        this parameter is deprecated as minimap2 is supplied in NextSV now
  -v, --version         show version number and exit

```


## Contact

If you have any questions/suggestions, please feel free to post it on the [Issue page](https://github.com/Nextomics/nextsv/issues). You may also email me at `fangli80@foxmail.com`. 

## Citation

### NextSV
Fang, L., Hu, J., Wang, D. et al. NextSV: a meta-caller for structural variants from low-coverage long-read sequencing data. BMC Bioinformatics 19, 180 (2018). https://doi.org/10.1186/s12859-018-2207-1

### Minimap2
Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100 (2018). https://doi.org/10.1093/bioinformatics/bty191

### Winnowmap
Jain, C., Rhie, A., Hansen, N.F., Koren, S., and Phillippy, A.M. (2022). Long-read mapping to repetitive reference sequences using Winnowmap2. Nat Methods 19, 705-710. https://doi.org/10.1038/s41592-022-01457-8

Jain, C., Rhie, A., Zhang, H., Chu, C., Walenz, B.P., Koren, S., and Phillippy, A.M. (2020). Weighted minimizer sampling improves long read mapping. Bioinformatics 36, i111-i118. https://doi.org/10.1093/bioinformatics/btaa435

### NGMLR and Sniffles
Sedlazeck, F.J., Rescheneder, P., Smolka, M. et al. Accurate detection of complex structural variations using single-molecule sequencing. Nat Methods 15, 461–468 (2018). https://doi.org/10.1038/s41592-018-0001-7

Smolka, M., Paulin, L.F., Grochowski, C.M., Mahmoud, M., Behera, S., Gandhi, M., Hong, K., Pehlivan, D., Scholz, S.W., Carvalho, C.M.B., et al. (2022). Comprehensive Structural Variant Detection: From Mosaic to Population-Level. bioRxiv. https://doi.org/10.1101/2022.04.04.487055 

### CuteSV
Jiang, T., Liu, Y., Jiang, Y. et al. Long-read-based human genomic structural variation detection with cuteSV. Genome Biol 21, 189 (2020). https://doi.org/10.1186/s13059-020-02107-y

Cao, S., Jiang, T., Liu, Y., Liu, S., and Wang, Y. (2022). Re-genotyping structural variants through an accurate force-calling method. bioRxiv. https://doi.org/10.1101/2022.04.04.487055 

### Samtools
Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

## Open source software tools included in this repository

To provide consistent output and better user experience, we included some open source software tools. 

**LongReadQC** https://github.com/fangli80/longreadqc

**Seqtk**
https://github.com/lh3/seqtk

**Minimap2**
https://github.com/lh3/minimap2

**NGMLR**
https://github.com/philres/ngmlr

**pigz**
https://github.com/madler/pigz

**Winnowmap**
https://github.com/marbl/Winnowmap
