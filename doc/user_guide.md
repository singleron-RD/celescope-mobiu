# User guide

## Requirements

- 64 bit Linux
- Minimum 32GB RAM to run [STAR](https://github.com/alexdobin/STAR) with human/mouse genome

## Installation

### Docker image

https://quay.io/repository/singleron-rd/celescope-mobiu?tab=tags

If you are unable to use Docker, you can install it using mamba/conda and pip as follows.

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
 --genomeDir /genome/rna/celescope_v2/hs \
 --kbDir /genome/kb/hs_ensembl99 \
 --thread 16 \
 --mod shell
```

`mapfile` Required. mapfile is a tab-delimited text file with four columns. Each line represents a pair of FASTQ files for a single sample.

| Column | Description |
|--------|-------------|
| 1st    | FASTQ file prefix (without `_R1`/`_R2` suffix). |
| 2nd    | Directory path to the FASTQ files. |
| 3rd    | Sample name, used as the prefix for all output files. |
| 4th    | Read end type: either `3p` for 3′-end sequencing or `5p` for 5′-end sequencing. |

Example

```tsv
3p_prefix   ./fastqs    test1   3p
5p_prefix   ./fastqs    test1   5p
```

```
$ ls ./fastqs

3p_prefix_001_R1.fastq.gz 3p_prefix_001_R2.fastq.gz
3p_prefix_002_R1.fastq.gz 3p_prefix_002_R2.fastq.gz
5p_prefix_001_R1.fastq.gz 5p_prefix_001_R2.fastq.gz
5p_prefix_002_R1.fastq.gz 5p_prefix_002_R2.fastq.gz
```

`--thread` It is recommended to use 16 threads. Using more than 20 threads is not advised because  [the mapping speed of STAR saturates at >~20 threads](https://github.com/singleron-RD/CeleScope/issues/197).

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

Start the analysis by running:
```
bash ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## Main output

- `outs/SJ.raw`  
  STARsolo raw splice junction matrix in [MARVEL](https://github.com/wenweixiong/MARVEL) input format. Chromosomes have been prefixed with "chr".

- `outs/raw`  
  Gene expression matrix file containing all barcodes (background + cells) from the barcode whitelist.

- `outs/filtered`  
  Gene expression matrix file containing only cell barcodes.

- `outs/transcript.filtered`  
  Transcript expression matrix containing only cell barcodes. Generated using [kb-python](https://github.com/pachterlab/kb_python).

- `outs/{sample}_Aligned.sortedByCoord.out.bam`  
  BAM file containing coordinate-sorted reads aligned to the genome. Read names starting with "5p:" indicate 5-prime reads, while all others are 3-prime reads.


## Seurat CreateSeuratObject
When using [seurat CreateSeuratObject](https://www.rdocumentation.org/packages/Seurat/versions/3.0.1/topics/CreateSeuratObject), the default `names.delim` is underscore . Since cell barcode is separated by underscore(for example, ATCGATCGA_ATCGATCGA_ATCGATCGA), using `names.delim = "_"` will incorrectly set `orig.ident` to the third segment of barcode. This problem can be avoided by setting names.delim to other characters, such as `names.delim="-"`
```
seurat.object = CreateSeuratObject(matrix, names.delim="-", project="sample_name") 
```

## Downstream analysis using [MARVEL](https://github.com/wenweixiong/MARVEL)
An example can be found at [MARVEL.ipynb](../jupyter/MARVEL.ipynb)

## [Change log](./CHANGELOG.md)
