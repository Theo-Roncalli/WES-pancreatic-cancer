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

# Checking if the install.sh script has been running before.
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

# What fraction of reads are mapped ? Almost 100%
# What fraction of pairs are properly mapped? Around 99.5%

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
echo -e "${BLUE}Running varscan somatic...${NC}";
varscan somatic ${mapping}/TCRBOA7-N-WEX-chr16.pileup ${mapping}/TCRBOA7-T-WEX-chr16.pileup ${variants}/TCRBOA7-T-WEX-chr16
echo -e "${GREEN}Done.${NC}\n";
# By using pileup files from a tumor-normal pair, varscan somatic gives some information about mutations.
# Two files are created: Variants.snp for single nucleotide polymorphism (SNP)
# and Variants.indel for insertions or deletions (indels).
# In these files, we can see the loss of heterozygosity (LOH), germline variants and somatic mutations.

# Checking if the script has been executed successfully.
if [ -f ${variants}/TCRBOA7-T-WEX-chr16.snp ] && [ -f ${variants}/TCRBOA7-T-WEX-chr16.indel ]
then
    echo -e "${GREEN}Success!";
	echo -e "You can find information about mutations by using the following files:"
    echo -e "${variants}/TCRBOA7-T-WEX-chr16.snp and ${variants}/TCRBOA7-T-WEX-chr16.indel${NC}"
    exit 0
else
    echo -e "\n${RED}The creation of files containing information about mutations has not been successully performed."
	echo -e "Please retry.${NC}"
    exit 1
fi
