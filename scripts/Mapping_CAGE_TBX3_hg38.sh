#!/bin/bash 

export PATH=/data/Shared/ruslan/hisat2-2.1.0:$PATH

reference='/data5/ruslan/human_heart_2/genome/ucsc' # path to the folder with genome .fasta file
ht_index='hg38_primary_chromosomes_ht2_index' # hisat index name
genome_fasta='hg38_primary_chromosomes.fa' # genome .fasta file itself, here I use only primary chromosomes
path_in='/home/shamil/human_muscles/final_fastq/first40' # path to final trimmed .fastq files
path_out='/home/shamil/human_muscles/mapped/hg38' # output directory

tr=20  # number of cores to use

#loop

for i in $path_in/*.fastq; do {

res_file=${i%.fastq*}.bwa.sai
res_file=${res_file/$path_in/}
res_file=$path_out$res_file
#echo $res_file
#} ; done

# BWA alignment

bwa aln -n 0.02 -o 1 -e -1 -i 5 -d 10 -l 32 -k 2 -m 2000000 -t $tr -M 3 -O 11 -E 4 -R 30 -q 0 $reference/$genome_fasta $i > $res_file;
bwa samse -n 100 $reference/$genome_fasta $res_file $i> ${res_file%.bwa.sai*}.sam;
samtools view --threads $tr -bSo ${res_file%.bwa.sai*}.bam ${res_file%.bwa.sai*}.sam;
samtools sort --threads $tr -m 3G ${res_file%.bwa.sai*}.bam -o ${res_file%.bwa.sai*}.sorted;
samtools index  -@ $tr ${res_file%.bwa.sai*}.sorted;
samtools view --threads $tr -b -F 4 ${res_file%.bwa.sai*}.sorted > ${res_file%.bwa.sai*}.mapped.bam;
samtools view --threads $tr -b -f 4 ${res_file%.bwa.sai*}.sorted > ${res_file%.bwa.sai*}.bwa_unmapped.bam;
bamToFastq -i ${res_file%.bwa.sai*}.bwa_unmapped.bam -fq ${res_file%.bwa.sai*}.bwa_unmapped.fastq;

# Hisat2 re-alignment of unmapped reads

hisat2 --threads $tr -x $reference/$ht_index -U ${res_file%.bwa.sai*}.bwa_unmapped.fastq  -S ${res_file%.bwa.sai*}.hisat2.sam;
samtools view --threads $tr -bSo ${res_file%.bwa.sai*}.hisat2.bam ${res_file%.bwa.sai*}.hisat2.sam;
samtools merge --threads $tr ${res_file%.bwa.sai*}.merged.bam ${res_file%.bwa.sai*}.hisat2.bam ${res_file%.bwa.sai*}.mapped.bam;
samtools sort --threads $tr -m 3G ${res_file%.bwa.sai*}.merged.bam -o ${res_file%.bwa.sai*}.merged.sorted.bam
samtools flagstat --threads $tr ${res_file%.bwa.sai*}.merged.sorted.bam > ${res_file%.bwa.sai*}.merged.sorted_flagstat.txt
samtools view --threads $tr -b -F 4 ${res_file%.bwa.sai*}.merged.sorted.bam > ${res_file%.bwa.sai*}.merged.sorted_mapped.bam;
samtools view --threads $tr -b -f 4 ${res_file%.bwa.sai*}.merged.sorted.bam > ${res_file%.bwa.sai*}.merged.sorted_unmapped.bam;

# CTSS.osc file (https://zenbu-wiki.gsc.riken.jp/zenbu/wiki/index.php/OSCtable)

# level1.py (http://genome.gsc.riken.jp/plessy-20150516/PromoterPipeline_20150516.tar.gz)

python2 /home/ruslan/Ruslan/ENC/ruslan/school/scripts/level1.py -o ${res_file%.bwa.sai*}.merged.sorted_mapped.osc ${res_file%.bwa.sai*}.merged.sorted_mapped.bam;

} ; done
