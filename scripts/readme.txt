here I attached:
1. A .docx file - "trimming_pipeline_CAGE.docx" with explanations how I do trimming; some scripts I took from Moirai pipeline (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-144), so I attached them in the archive file
2. scripts for trimming and alignment - "scripts.zip"

Mapping_CAGE_TBX3_hg38.sh - use for alignment of trimmed reads to reference genome assembly. There are several variables that need to be specified, like 'path_in', 'path_out' etc.

If everything goes well you can get .osc files - tables with CAGE signal ( CTSS ) across the genome. These files can be used for promoters or enhancers calling:
1. promoters (http://genome.gsc.riken.jp/plessy-20150516/PromoterPipeline_20150516.tar.gz): python2 level2.py -o /path/to/output/promoters.osc /path/to/ctss/*.osc
2. enhancers: https://github.com/anderssonrobin/enhancers ; https://www.nature.com/articles/nature12787 ;  I will make :)

simple way for promoter to gene connection was in my R code I shared previously - "tbx3_annotation_DE.R"

Motif analysis I usually do with MEME-Suite:
https://meme-suite.org/meme/

For example AME https://meme-suite.org/meme/tools/ame:
+300 TSS -100 sequences in fasta format for DE CAGE peaks vs non-DE CAGE peaks as background with vertebrate motif collection

Let me know if you have any questions or something doesn't work!

=======================================
There are several pipelines for CAGE analysis exist:

CAGEr (http://bioconductor.org/packages/release/bioc/html/CAGEr.html) - useful for TSS shifts
CAGEfightR (http://bioconductor.org/packages/release/bioc/html/CAGEfightR.html) - maybe can predict enhancers, I haven't try yet
PromoterPipeline  (http://genome.gsc.riken.jp/plessy-20150516/PromoterPipeline_20150516.tar.gz) - I use for promoters calling, very simple python script
PDI ( https://github.com/hkawaji/dpi1/ ) - classic FANTOM5 protocol for promoters

Regards,
Ruslan

