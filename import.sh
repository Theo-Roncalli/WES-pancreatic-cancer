#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads

# Url parameters

reads_url=https://drive.google.com/drive/folders/1R4DEQ56jPY3wior5Mow359oPplDrZKD-/patient7.exome.tar.gz

# Step 1: Download reads

mkdir -p ${reads}
echo "Downloading reads..."
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" \
    -O ${reads}/patient7.tar.gz && rm -rf /tmp/cookies.txt # The options in the command wget allow to download a large file on Google Drive.

