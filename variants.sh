#!/usr/bin/env bash

# Directory parameters

trimming=Data/Trimming
index=Data/Index
mapping=Data/Mapping
variants=Data/Variants

# Performance parameters

nb_cpus_mapping=7

# Color parameters

RED='\033[1;31m'
GREEN='\033[1;32m'
BLUE='\033[1;34m'
NC='\033[0m'

if [ ! -d ${trimming} ] || [ ! -d ${index} ]
then
    echo -e "${RED}Import data must be done before running the variants.sh script.";
    echo -e "Please run install.sh script before with the command: \033[3;31mbash install.sh${NC}";
    exit 1
fi

# Step 1: Mapping

#### Mapping (BWA) ####
mkdir ${mapping} -p
for read1_file in ${trimming}/*1P.fastq
do
	paired_file_with_path=${read1_file%_1P.fastq};
	paired_file_without_path=${paired_file_with_path#${trimming}/};
	echo -e "${BLUE}Downloading SAM file with ${paired_file_without_path}...${NC}";
    bwa mem -M -t ${nb_cpus_mapping} -A 2 -E 1 \
        ${index}/chr16.fa.gz \
        ${paired_file_with_path}_1P.fastq \
        ${paired_file_with_path}_2P.fastq \
        -o ${mapping}/${paired_file_without_path}.sam;
    echo -e "${GREEN}Done.${NC}\n";
done

# Options used:
# -t INT        Number of threads 
# -A INT        Matching score. 
# -E INT        Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap).
# -M            Mark shorter split hits as secondary (for Picard compatibility)
# -O INT        Gap open penalty.

# Step 2: Processing

for samfile in ${mapping}/*.sam
do
    bamfile=${samfile%.sam}.bam;
    sorted_bamfile=${samfile%.sam}.sorted.bam;
    echo -e "${BLUE}Converting ${samfile} into ${bamfile}...${NC}";
    samtools view -S -b ${samfile} > ${bamfile};            # converting SAM to BAM
    echo -e "\n${BLUE}Checking Mapping Stats...${NC}";
    samtools flagstat ${bamfile};                           # checking statistics
    echo -e "\n${BLUE}Sorting ${bamfile}...${NC}";
    samtools sort ${bamfile} > ${sorted_bamfile};           # sorting BAM file
    echo -e "\n${BLUE}Indexing ${sorted_bamfile}...${NC}";
    samtools index ${sorted_bamfile};                       # indexing BAM file
    echo -e "${GREEN}Done.${NC}\n";
done

# What fraction of reads were mapped ? Around 99.98%
# What fraction of pairs were properly mapped? Around 0%

gunzip ${index}/chr16.fa.gz
for sorted_bamfile in ${mapping}/*.sorted.bam
do
    pileup_file=${sorted_bamfile%.sorted.bam}.pileup;
    echo -e "${BLUE}Creating ${pileup_file}...${NC}";
    samtools mpileup -B -A -f ${index}/chr16.fa ${sorted_bamfile} > ${pileup_file}; # converting to pileup
    echo -e "${GREEN}Done.${NC}\n";
done

# Step 3: Calling somatic variants with Varscan

# sudo apt install varscan
mkdir ${variants} -p
# varscan somatic ${mapping}/*-N-*.pileup ${mapping}/*-T-*.pileup ${variants}
varscan somatic ${mapping}/*-N-WEX-chr16.pileup ${mapping}/*-T-WEX-chr16.pileup ${variants}

exit