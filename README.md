# Trimming

### Merging fastq files
Use "merge_fastq.sh" to combine fastq files by lanes.

### Renaming
Rename files to add sequence information via "new_names.tsv" file.
`awk -F'\t' '{gsub(/[[:space:]]+$/, "", $2); gsub(/^'"'"'|'"'"'/, "", $2); gsub(/[$'\r']$/, "", $2); printf "mv \"%s\" \"%s\"\n", $1, $2}' new_names.tsv > rename_commands.sh`

### Count the reads
`for fastq in *.fastq; do echo "$(basename $fastq) $(expr $(wc -l < $fastq) / 4)" >> number_of_***_reads.txt ; done`
