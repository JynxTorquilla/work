
# TBX3

setwd("E:/work/tbx3")

osc_tab <- read.table("./tbx3_level2.osc", header = T)

library(scales)
library(factoextra)
library(RColorBrewer)
library(edgeR)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(GenomicFeatures)
library(ChIPseeker)
library(AnnotationDbi)
library(org.Hs.eg.db)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

PromotersTSS <- GRanges(seqnames = osc_tab$chrom, ranges = IRanges(start =osc_tab$start.0base,end = osc_tab$end),strand = osc_tab$strand, name=osc_tab$id, tss = osc_tab$pos)
PromotersTSS

annotatePeak_hg38 <- annotatePeak(PromotersTSS, tssRegion=c(-3000, 3000), TxDb=txdb, sameStrand = T)
annotatePeak_hg38

plotAnnoPie(annotatePeak_hg38)

as.data.frame(annotatePeak_hg38@anno) -> ann_cage_hs
ann_cage_hs$GeneName <- egSYMBOL$symbol[match(ann_cage_hs$geneId, egSYMBOL$gene_id)]

all(osc_tab$id == ann_cage_hs$name)

exp_tab <- osc_tab[,grep("norm.", colnames(osc_tab))]
colnames(exp_tab)

group <- as.factor(c(rep("E", 3), rep("T", 2)))
group

c <- DGEList(counts = exp_tab, genes = ann_cage_hs, group = group)
design <- model.matrix(~0+group)
c <- calcNormFactors(c, method = "TMM")
c <- estimateDisp(c)
c <- estimateCommonDisp(c)
c <- estimateTagwiseDisp(c)
norm_counts.table <- t(t(c$pseudo.counts)*(c$samples$norm.factors))

res.pca <-prcomp(t(norm_counts.table), scale = T)
p5<-fviz_pca_ind(res.pca,
                 col.ind = group, # Color by the quality of representation
                 #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 geom = c("point"),
                 title = "tmm",
                 show.legend = T,
                 repel = T, addEllipses = T
)
plot(p5)

colors <- brewer.pal(12,"Set3")[c(4,5,6,7,8)]
points <- c(15,16,17,18,19)
md_data<-plotMDS(c, col=alpha(colors[group],0.8), pch=points[group], main= '', cex=1.5, ylab="", xlab="") # , top = dim(exp_tab)[1]

et <- exactTest(c)
et$comparison
topTags(et, n = nrow(exp_tab)) -> de_tab_1
de_tab_1 <- de_tab_1$table

dim(de_tab_1[de_tab_1$FDR < 0.05, ])

de_tab_1->final1

final1$logPval<-abs(log10(final1$PValue))
final1$Colour<-"black"
final1$Colour[final1$FDR<0.05 & final1$logFC> 1]<-brewer.pal(12, "Paired")[10] # orange 8 6
final1$Colour[final1$FDR<0.05 & final1$logFC< -1]<-brewer.pal(12, "Paired")[12] # blue 2 4

head(final1)

n_genes <- 20
final1$labels <- NA
final1 <- final1[order(final1$logFC * final1$logPval, decreasing = T), ]
final1$labels[!is.na(final1$GeneName) & final1$FDR<0.05 & final1$logFC> 1 & final1$GeneName!="" & !duplicated(final1$GeneName)][1:n_genes] <- final1$GeneName[!is.na(final1$GeneName) & final1$FDR<0.05 & final1$logFC> 1  &final1$GeneName!="" & !duplicated(final1$GeneName)][1:n_genes]
final1 <- final1[order(final1$logFC * final1$logPval, decreasing = F), ]
final1$labels[!is.na(final1$GeneName) & final1$FDR<0.05 & final1$logFC < -1 &  final1$GeneName!="" &!duplicated(final1$GeneName)][1:n_genes] <- final1$GeneName[!is.na(final1$GeneName) & final1$FDR<0.05 & final1$logFC < -1 &final1$GeneName!="" & !duplicated(final1$GeneName)][1:n_genes]
length(unique(final1$GeneName[final1$FDR < 0.05 & final1$logFC > 0]))
length(unique(final1$GeneName[final1$FDR < 0.05 & final1$logFC < 0]))

library(ggrepel)
colnames(final1)
colnames(final1)[duplicated(colnames(final1))]
ma <- ggplot(data = final1, aes(y = logPval , x = logFC)) +   geom_point(colour=final1$Colour) + theme_light()   + geom_vline(xintercept=c(-1,1), linetype="dashed", color = "gold")  +  geom_label_repel(aes(label=labels),  data=final1,max.overlaps = 5000, label.size = NA,  fill=NA) + ggtitle(paste("down: ", length(unique(final1$GeneName[final1$FDR < 0.05 & final1$logFC < 0])), "; up: ", length(unique(final1$GeneName[final1$FDR < 0.05 & final1$logFC > 0])), sep=""))
ma

final1_counts <- merge(final1, osc_tab[,c(1,7:16)], by.x="name", by.y="id", all.x=T)

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16192m"))
library(xlsx)
final1_counts <- final1_counts[order(final1_counts$logFC * final1_counts$logPval, decreasing = T), ]
write.xlsx2(x= final1_counts, file = "TBX3 CAGE annotation, DE, counts.xlsx", col.names = T, row.names = F)


gene_list <- final1_counts$logFC[final1_counts$PValue < 0.05]
names(gene_list) <- final1_counts$geneId[final1_counts$PValue < 0.05]

gene_list <- gene_list[!is.na(names(gene_list))]
gene_list = sort(gene_list, decreasing = TRUE)
length(gene_list)

library(clusterProfiler)
library(viridis)
library(org.Hs.eg.db)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH",
             eps = 0)


go_plot_1 <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) +scale_color_viridis(option = "A", begin= 0.5, end = 0.9)# + scale_color_gradient(high="#FF7F00", low= "#1F78B4")
go_plot_1


kegg_organism = "hsa"
kk2_a <- gseKEGG(geneList     = gene_list,
                 organism     = kegg_organism,
                 #nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 1000,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 eps = 0)

kegg_plot_1 <- dotplot(kk2_a, showCategory=10, split=".sign") + facet_grid(.~.sign) +scale_color_viridis(option = "A", begin= 0.5, end = 0.9)# + scale_color_gradient(high="#FF7F00", low= "#1F78B4")
kegg_plot_1

pdf("./tbx3 CAGE figures.pdf", width = 7, height = 7)
plotAnnoPie(annotatePeak_hg38, main = "CAGE peaks distribution") 
p5
ma
go_plot_1
kegg_plot_1
dev.off()



