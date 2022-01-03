#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads
trimming=Data/Trimming
genome=Data/Genome
index=Data/Index
mapping=Data/Mapping
counts=Data/Counts
figures_reads=Figures/Reads
figures_trimming=Figures/Trimming

# Performance parameters

nb_cpus_indexing=7
nb_cpus_mapping=7
nb_cpus_counting=7

# Step 1: Quality control + Reads cleaning

mkdir -p ${figures_reads}
echo "Creation of the fastqc files on raw reads."
fastqc -o ${figures_reads} -f fastq ${reads}/*.fastq -q

# Trimming procedure (Elimination of low quality sequences at the end of reads)
# conda install -c bioconda trimmomatic

mkdir -p ${trimming}

for read1_file in ${reads}/*_r1F.fastq
do
	paired_file_with_path=${read1_file%_r1F.fastq};
	paired_file_without_path=${paired_file_with_path#${reads}/};
	echo "Trimming ${paired_file_without_path%.sampled}...";
    echo ${paired_file_with_path}_r1F.fastq
	trimmomatic PE ${paired_file_with_path}_r1F.fastq ${paired_file_with_path}_r2F.fastq -baseout ${trimming}/${paired_file_without_path}.fastq LEADING:20 TRAILING:20 MINLEN:50 -quiet
	echo "Done."
done

# Removal of the bases from the extremity with a quality lower than 20. If the final read is smaller than 50, it is discarded. file with U => discard. file with P => no discard.
# Remark: files 1U and 2U returns a small number of sequences (around 10 000) while files 1P and 2P returns a large number of sequences (a little smaller than the reads without cleaning)

mkdir -p ${figures_trimming}
echo "Creation of the fastqc files on trimmed reads."
fastqc -o ${figures_trimming} -f fastq ${trimming}/*.fastq -q

