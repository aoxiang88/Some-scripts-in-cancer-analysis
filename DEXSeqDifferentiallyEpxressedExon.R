# for computing the differential usage of exon
.libPaths(c(.libPaths(), "/home/aoxiang/lib/R", "/home/aoxiang/_bio_Tools/R-3.4.3/library"))
library("DEXSeq")
library("BiocParallel")
library("hwriter")
source("/home/aoxiang/_bio_Tools/subread_to_DEXSeq/load_SubreadOutput.R")
setwd("/home/aoxiang/data/bladder-cancer/clean_fq/bam_files_by_STAR/")
clinical_list <- read.table("clinical_characteristics_RNAseq.list", header=TRUE)
clinical_frame <- data.frame(row.names=clinical_list$CaseID, condition=clinical_list$State, age=clinical_list$Age, sex=clinical_list$Sex, grade=clinical_list$Grade, position=clinical_list$Superficial.Invasive, primary=clinical_list$Primary.Relapsed)
dxd.fc <- DEXSeqDataSetFromFeatureCounts("readsCounts_4_exons_58.txt", flattenedfile="bladder_all_merged_4_featureCounts.gtf", sampleData=clinical_frame)
dxd <- dxd.fc[rowSums(featureCounts(dxd.fc))>0,] # filter out the exons that were not expressed at all 
formulaFullModel <- ~ sample + exon + age:exon + sex:exon + grade:exon + position:exon + primary:exon + condition:exon
formulaReducedModel <- ~ sample + exon + age:exon + sex:exon + grade:exon + position:exon + primary:exon 
BPPARAM <- MulticoreParam(worker=5)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
#pdf("exon_qc_dispersions.pdf",10,10)
#plotDispEsts(dxd)
#dev.off()
dxd = testForDEU( dxd, BPPARAM=BPPARAM, reducedModel = formulaReducedModel, fullModel = formulaFullModel)
dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
dxr2 = DEXSeqResults(dxd)
#DEXSeqHTML( dxr2, FDR=0.01, color=c("#FF000080", "#0000FF80") )
#pdf("exon_qc_dispersions.pdf",10,10)
#plotDispEsts(dxd)
#dev.off()
#pdf("exon_qc_MA_plot.pdf",10,10)
#plotMA(dxr2, alpha=0.01, cex=0.8)
#dev.off()
save.image("exon_differential_usage.RData")
