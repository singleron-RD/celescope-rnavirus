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
wget https://raw.githubusercontent.com/singleron-RD/celescope-mobiu/refs/heads/master/conda_pkgs.txt?token=GHSAT0AAAAAACQ7BHHCKNDZNVI4I6WWVAQ4Z74XZBA
```

Then, start creating the conda environment and install the conda package.It is recommended to use [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (which is a faster replacement for Conda) to install conda packages.
The following command will create a conda environment named `celescope` and install the dependency packages contained in conda_pkgs.txt.
```
mamba create -n mobiu -y --file conda_pkgs.txt
```

### Install celescope-mobiu

Make sure you have activated the conda environment before running `pip install`.
```
mamba activate mobiu
pip install celescope-mobiu
```

## Usage

1. Create a genomeDir.

### Homo sapiens

```
mkdir hs_ensembl_99
cd hs_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz

conda activate celescope
celescope utils mkgtf Homo_sapiens.GRCh38.99.gtf Homo_sapiens.GRCh38.99.filtered.gtf
celescope rna mkref \
 --genome_name Homo_sapiens_ensembl_99_filtered \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.99.filtered.gtf \
 --mt_gene_list mt_gene_list.txt
```

### Mus musculus

```
mkdir mmu_ensembl_99
cd mmu_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz 
gunzip Mus_musculus.GRCm38.99.gtf.gz

conda activate celescope
celescope utils mkgtf Mus_musculus.GRCm38.99.gtf Mus_musculus.GRCm38.99.filtered.gtf

celescope rna mkref \
 --genome_name Mus_musculus_ensembl_99_filtered \
 --fasta Mus_musculus.GRCm38.dna.primary_assembly.fa \
 --gtf Mus_musculus.GRCm38.99.filtered.gtf \
 --mt_gene_list mt_gene_list.txt
```

2. Create a [kb reference](https://github.com/pachterlab/kb_python?tab=readme-ov-file#kb-ref-generate-a-pseudoalignment-index) directory.

```
kb ref \
 --workflow=standard \
 -i index.idx -g t2g.txt -f1 cdna.fa \
 Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.99.filtered.gtf
```

3. Generate scripts for each sample.

Under your working directory, write a shell script `run.sh` as
```
multi_mobiu \
 --mapfile mapfile \
 --chemistry mobiu-1 \
 --genomeDir /genome/rna/celescope_v2/hs \
 --kbDir /genome/kb/hs_ensembl99 \
```

**mapfile Format**
```tsv
fastq_prefix_5p fastq_5p_folder sample_name 5p
fastq_prefix_3p fastq_3p_folder sample_name 3p
```

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## [Change log](./CHANGELOG.md)




 
