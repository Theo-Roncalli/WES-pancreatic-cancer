#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads

# Url parameters

reads_url="'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe'"

# Step 1: Download reads

mkdir -p ${reads}
echo "Downloading reads..."
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" \
    -O ${reads}/patient7.tar.gz && rm -rf /tmp/cookies.txt # The options in the command wget allow to download a large file on Google Drive.
echo "Unarchiving reads..."
tar -zxf ${reads}/patient7.tar.gz -C ${reads}
mv ${reads}/patient7.exome/* ${reads} && rm -r ${reads}/patient7.exome && rm ${reads}/patient7.tar.gz # Cleaning 'reads' repository
gunzip -q ${reads}/*.fastq.gz
