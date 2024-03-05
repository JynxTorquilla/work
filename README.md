# Trimming

### Merging fastq files
Use "merge_fastq.sh" to combine fastq files by lanes.

### Renaming
Rename files to add sequence information via "new_names.tsv" file.
awk -F'\t' '{gsub(/[[:space:]]+$/, "", $2); gsub(/^'"'"'|'"'"'/, "", $2); gsub(/[$'\r']$/, "", $2); printf "mv \"%s\" \"%s\"\n", $1, $2}' new_names.tsv > rename_commands.sh

