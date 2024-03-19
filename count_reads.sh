#!/bin/bash

# Create a new text file to store the results
echo "Filename Reads" > number_of_reads.txt

# Loop through all fastq.gz files in the current directory
for fastq in *.fastq.gz; do
  # Count the number of reads in the file
  reads=$(zcat $fastq | wc -l | awk '{print $1/4}')
  # Write the filename and number of reads to the text file
  echo "$fastq $reads" >> number_of_reads.txt
done