#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads
trimming=Data/Trimming
genome=Data/Genome
index=Data/Index
figures_reads=Figures/Reads
figures_trimming=Figures/Trimming

# Url parameters

genome_url=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz

# Color parameters

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Step 1: Download reads

mkdir -p ${reads}
echo "Downloading reads..."
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" \
    -O ${reads}/patient7.tar.gz && rm -rf /tmp/cookies.txt # The options in the command wget allow to download a large file on Google Drive.
echo "Unarchiving reads..."
tar -zxf ${reads}/patient7.tar.gz -C ${reads}
mv ${reads}/patient7.exome/* ${reads} && rm -r ${reads}/patient7.exome && rm ${reads}/patient7.tar.gz # Cleaning 'reads' repository
gunzip -q ${reads}/*.fastq.gz

# Print read numbers and file sizes for each sample
echo "--------Number of reads and size of each file--------"
for file in ${reads}/*.fastq
do
    # Print filename, number of reads and size of a specific file
	grep ^+$ ${file} | echo "$(wc -l) $(du -h ${file})" | awk '{print $3 " Reads:" $1 " Size:" $2}';
done

# Check if read numbers are the same between paired files
# Exit if not the case
echo -e "\n"
for read1_file in ${reads}/*_r1F.fastq
do
	paired_file_with_path=${read1_file%_r1F.fastq};
    if [ $(grep ^+$ ${paired_file_with_path}_r1F.fastq | wc -l) == $(grep ^+$ ${paired_file_with_path}_r2F.fastq | wc -l) ]
    then
        echo -e "${GREEN}Number of reads between the paired files ${paired_file_with_path}_r1F.fastq and ${paired_file_with_path}_r2F.fastq are the same.${NC}"
    else
        echo -e "${RED}Error: number of reads between the paired files ${paired_file_with_path}_r1F.fastq and ${paired_file_with_path}_r2F.fastq are not the same.${NC}"
        exit 1;
    fi
done
echo -e "\n"

# Step 2: Quality control + Reads cleaning

mkdir -p ${figures_reads}
echo "Creation of the fastqc files on raw reads."
fastqc -o ${figures_reads} -f fastq ${reads}/*.fastq -q
echo "Done."

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

# Step 3: Download reference genome

mkdir ${index} -p
echo "Downloading genome..."
wget ${genome_url} -P ${index} -q

# Step 4: Creation of the index

#### Index (BWA) ####
mkdir ${index} -p
bwa index -a bwtsw ${index}/chr16.fa.gz
# mv ${genome}/*.fa.gz.+([A-Za-z]) ${index}/
