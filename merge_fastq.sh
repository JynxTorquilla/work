#!/bin/bash

for file in *__L00*.fastq.gz; do
    prefix="${file%%__L00*}"  # Extract the part before "__L00"
    suffix="${file#*__L00}"   # Extract the part after "__L00"
    suffix="${suffix%.fastq.gz}"  # Remove ".fastq.gz" from the end
    cat "${prefix}"*__L00*.fastq.gz > ./merged_fastq/"${prefix}.fastq.gz" 
done