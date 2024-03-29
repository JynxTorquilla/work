#!/bin/bash

mkdir merged_fastq_2

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > ./merged_fastq_2/"$i"_L001_R1_001_merged.fastq.gz

done;
