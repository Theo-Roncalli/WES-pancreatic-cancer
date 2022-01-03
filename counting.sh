#!/usr/bin/env bash

# Directory parameters

# reads=Data/Reads
trimming=Data/Trimming
index=Data/Index
# mapping=Data/Mapping
# counts=Data/Counts
mapping=Data/Mapping

# Performance parameters

nb_cpus_mapping=7

# Color parameters

RED='\033[1;31m'
GREEN='\033[1;32m'
BLUE='\033[1;34m'
NC='\033[0m'

# Step 1: Mapping

#### Mapping (BWA) ####
mkdir ${mapping} -p
for read1_file in ${trimming}/*1P.fastq
do
	paired_file_with_path=${read1_file%_1P.fastq};
	paired_file_without_path=${paired_file_with_path#${trimming}/};
	echo -e "${BLUE}Downloading SAM file with ${paired_file_without_path}...${NC}\n";
    bwa mem -M -t ${nb_cpus_mapping} -A 2 -E 1 \
        ${index}/chr16.fa.gz \
        ${paired_file_with_path}_1P.fastq \
        ${paired_file_with_path}_1P.fastq \
        -o ${mapping}/${paired_file_without_path}.sam;
    echo -e "${GREEN}Done.${NC}\n";
done

# Options used:
# -t INT        Number of threads 
# -A INT        Matching score. 
# -E INT        Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap).
# -M            Mark shorter split hits as secondary (for Picard compatibility)

# Step 2: Processing

for samfile in ${mapping}/*.sam
do
    bamfile=${samfile%.sam}.bam;
    sorted_bamfile=${samfile%.sam}.sorted.bam;
    echo -e "${BLUE}Converting ${samfile} into ${bamfile}...${NC}";
    samtools view -S -b ${samfile} > ${bamfile};            # converting SAM to BAM
    echo -e "\n${BLUE}Checking Mapping Stats...${NC}";
    samtools flagstat ${bamfile};                           # checking statistiques
    echo -e "\n${BLUE}Sorting ${bamfile}...${NC}";
    samtools sort ${bamfile} > ${sorted_bamfile};           # sorting BAM file
    echo -e "\n${BLUE}Indexing ${sorted_bamfile}...${NC}";
    samtools index ${sorted_bamfile};                       # indexing BAM file
    echo -e "${GREEN}Done.${NC}\n";
done
