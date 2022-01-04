# New Generation Sequencing (NGS)

This work focuses on the study of [Becnel et al. (2016)](https://www.nature.com/articles/sdata201610) who share cancer genomic data. In our work, we are interesting into mutations that could induce pancreatic cancer. Here, we work on data from patient ID TCRBOA7. In order to reduce time computation, we only used data on chromosome 16.

## Dependencies

The pipeline runs on bash.
Some package are required for launching some commands such as fastqc, trimmomatic, bwa and varscan.

```bash
sudo apt install fastqc                 # For using fastqc
conda install -c bioconda trimmomatic   # For using trimmomatic
sudo apt install bwa                    # For using the bwa library
sudo apt install varscan                # For using varscan somatic
```

## Hardware requirements

A machine with at least 8 GB of **FREE** RAM (to create the index and the mapping on the chromosome 16 of the reference genome).

## Executing The Pipeline

The pipeline is used for detecting variants on whole exome sequencing (WES) data using paired-end files: the first two for tumor tissues and the other two for adjacent normal tissues. Here, it is used on a patient suffering from pancreatic cancer. Two files are generated at the end of the exome sequencing pipeline: a file with the single nucleotide polymorphisms (SNP) and a file with the insertions/deletions (indels). These files are available in the repository _Data/Variants_. The steps for generating these files are the followings.

1. Clone the Github repository to your machine
```bash
git clone https://github.com/Theo-Roncalli/WES-pancreatic-cancer.git
cd WES-pancreatic-cancer
```

2. Importation of reads and reference genome
```bash
bash install.sh
```

3. Creation of the variant files which contains the SNP and indels.
```bash
bash variants.sh
```

## Cleaning Repository

For cleaning the repository (i.e. delete _Data_ and _Figures_ folders), please type:
```bash
bash clean.sh
```