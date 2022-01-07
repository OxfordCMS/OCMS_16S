#!/bin/bash
input="download_files.tsv"
while IFS= read -r line
do
  wget $line
done < $input
