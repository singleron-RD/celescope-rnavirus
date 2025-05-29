# User guide

## Hardware/Software requirements

- 64 bit Linux
- Minimum 32GB RAM to run [STAR](https://github.com/alexdobin/STAR) with human/mouse genome
- [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or [conda](https://anaconda.org/anaconda/conda)

## Installation

### Create conda environment and install conda packages. 
First, you need to get the txt file from the github repository containing the name of the conda package you need to install. You can download it directly from the github repository interface, or use a download link.
The following command will download the conda_pkgs.txt required for the latest version.
```
wget https://raw.githubusercontent.com/singleron-RD/celescope-rnavirus/refs/heads/master/conda_pkgs.txt
```

Then, start creating the conda environment and install the conda package.It is recommended to use [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (which is a faster replacement for Conda) to install conda packages.
The following command will create a conda environment named `celescope` and install the dependency packages contained in conda_pkgs.txt.
```
mamba create -n rnavirus -y --file conda_pkgs.txt
```

### Install celescope-rnavirus

Make sure you have activated the conda environment before running `pip install`.
```
mamba activate rnavirus
pip install celescope-rnavirus
```

## Usage

1. Download a kb reference.

https://github.com/pachterlab/LSCHWCP_2023/tree/main/precomputed_refs
https://github.com/pachterlab/LSCHWCP_2023/tree/main/PalmDB

2. Generate scripts for each sample.

Under your working directory, write a shell script `run.sh` as
```
multi_rnavirus \
 --mapfile mapfile \
 --kbDir /genome/kb/human_mask_ref \
 --mod shell
```

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## [Change log](./CHANGELOG.md)




 
