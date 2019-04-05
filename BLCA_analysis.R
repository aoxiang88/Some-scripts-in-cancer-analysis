# =====================================================================
# For compatibility with Rscript.exe: 
# =====================================================================
if(length(.libPaths()) == 1){
  # We're in Rscript.exe
  possible_lib_paths <- file.path(Sys.getenv(c('USERPROFILE','R_USER')),
                                  "R","win-library",
                                  paste(R.version$major,
                                        substr(R.version$minor,1,1),
                                        sep='.'))
  indx <- which(file.exists(possible_lib_paths))
  if(length(indx)){
    .libPaths(possible_lib_paths[indx[1]])
  }
  # CLEAN UP
  rm(indx,possible_lib_paths)
}
# =====================================================================
#source("https://bioconductor.org/biocLite.R") 
#biocLite("DESeq2")
#biocLite("tximport")
library(edgeR)
library(DESeq2)
library(apeglm)
library(tximport)
library(pheatmap)
setwd("/home/aoxiang/data/bladder-cancer/rsem_results/")
clinic <- read.table("../clean_fq/clinical_characteristics.txt", header = TRUE)
designMat <- data.frame(condition = clinic$State, patient = clinic$Patient, sample = clinic$CaseID, age = as.factor(clinic$Year), sex = clinic$Sex)
rownames(designMat) <- clinic$CaseID

# construct data from tumor samples that do not have paired normal
normalIndex <- which(designMat$condition == "Normal")
singleTumorIndex <- which(! designMat$patient %in% designMat$patient[normalIndex])
pairedNormalSample <- unlist(lapply(designMat$patient[singleTumorIndex], FUN = function(x){return(paste0(as.character(x), "-N"))}))
pairedNormalPatient <- designMat$patient[singleTumorIndex]
pairedNormalCondition <- rep("Normal", length(singleTumorIndex))
pairedNormalAge <- designMat$age[singleTumorIndex]
pairedNormalSex <- designMat$sex[singleTumorIndex]
pairedDesignMat <- data.frame(condition = pairedNormalCondition, patient = pairedNormalPatient, sample = pairedNormalSample, age = pairedNormalAge, sex = pairedNormalSex)
designMat <- rbind(designMat, pairedDesignMat)
rownames(designMat)[(nrow(clinic)+1): nrow(designMat)] <- pairedNormalSample

##-------- DEG analysis -------- ##
files <- file.path("./expression_data", paste0(clinic$CaseID, ".genes.results", sep=""))
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
colnames(txi$counts) <- clinic$CaseID

# edgeR
cts <- txi$counts
normMat <- txi$length
normMat[normMat == 0] <- 1
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y <- scaleOffset(y, t(t(log(normMat)) + o))
keep <- filterByExpr(y)
y <- y[keep, ]
pdf("./image/edgeR_MDS.pdf")
plotMDS(y)
dev.off()

# DESeq2
txi$length[txi$length==0] <- 1
genesNormalTPM <- read.table("./expression_data/genes.normal.TPM", sep = "\t", header = TRUE, row.names = 1)
genesTumorTPM <- read.table("./expression_data/genes.tumor.TPM", sep = "\t", header = TRUE, row.names = 1)
genesNormalMedian <- apply(genesNormalTPM, MARGIN = 1, median)
genesTumorMedian <- apply(genesTumorTPM, MARGIN = 1, median)
geneSelect <- which(genesNormalMedian > 1 | genesTumorMedian > 1) # filter out genes that were not expressed(median TPM < 1) under both conditions
txiFiltered <- list(abundance = txi$abundance[geneSelect, ], counts  = txi$counts[geneSelect, ], length = txi$length[geneSelect, ], countsFromAbundance = txi$countsFromAbundance)

normalCountMedian <- rowMedians(txiFiltered$counts[, normalIndex])
normalAbundMedian <- rowMedians(txiFiltered$abundance[, normalIndex])
normalLenMedian <- rowMedians(txiFiltered$length[, normalIndex])
normalCountMat <- rep(normalCountMedian, length(singleTumorIndex))
normalCountMat <- matrix(normalCountMat, length(normalCountMedian), length(singleTumorIndex))
colnames(normalCountMat) <- pairedNormalSample
normalAbundMat <- rep(normalAbundMedian, length(singleTumorIndex))
normalAbundMat <- matrix(normalAbundMat, length(normalAbundMedian), length(singleTumorIndex))
colnames(normalAbundMat) <- pairedNormalSample
normalLenMat <- rep(normalLenMedian, length(singleTumorIndex))
normalLenMat <- matrix(normalLenMat, length(normalLenMedian), length(singleTumorIndex))
colnames(normalLenMat) <- pairedNormalSample
txiFiltered$abundance <- cbind(txiFiltered$abundance, normalAbundMat)
txiFiltered$counts <- cbind(txiFiltered$counts, normalCountMat)
txiFiltered$length <- cbind(txiFiltered$length, normalLenMat)

dds <- DESeqDataSetFromTximport(txiFiltered, designMat, ~patient + condition)
#select <- rowSums(counts(dds)) >=10
#dds <- dds[select,]
dds$condition <- relevel(dds$condition, ref = "Normal")
dds <- DESeq(dds)
dds <- dds[which(mcols(dds)$betaConv),]
#resLFC <- lfcShrink(dds, coef = "condition_Tumor_vs_Normal", type = "apeglm")
res <- results(dds)
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), file = "./results/DEG.00.paired_results_by_FDR.csv")
select <- which(abs(resOrdered$log2FoldChange) > 2 & res$padj < 0.01)
vsd <- vst(dds, blind = FALSE)
ntd <- normTransform(dds)
rownames(anno) <- clinic$CaseID
anno <- data.frame(state = clinic$State)
pheatmap(assay(ntd)[select,1:58], filename = "./results/DEG.Paired_cluster.pdf", show_rownames = FALSE, annotation_col = anno, fontsize_col = 6, cutree_cols = 3)

#DET analysis
filesTx <- file.path("./expression_data", paste0(clinic$CaseID, ".isoforms.results", sep = ""))
txiTx <- tximport(filesTx, type = "rsem", txIn = TRUE, txOut = TRUE)
colnames(txiTx$counts) <- clinic$CaseID
txiTx$length[txiTx$length==0] <- 1
txNormalTPM <- read.table("./expression_data/tx.normal.TPM", header = TRUE, sep = "\t", row.names = 1)
txTumorTPM <- read.table("./expression_data/tx.tumor.TPM", header = TRUE, sep = "\t", row.names = 1)
txNormalMedian <- apply(txNormalTPM, MARGIN=1, median)
txTumorMedian <- apply(txTumorTPM, MARGIN=1, median)
txSelect <- which(txNormalMedian > 1 | txTumorMedian > 1)
txiTxFiltered <- list(abundance = txiTx$abundance[txSelect, ], counts = txiTx$counts[txSelect, ], length = txiTx$length[txSelect, ], countsFromAbundance = txiTx$countsFromAbundance)
# Recover
tx_normalCountMedian <- rowMedians(txiTxFiltered$counts[, normalIndex])
tx_normalAbundMedian <- rowMedians(txiTxFiltered$abundance[, normalIndex])
tx_normalLenMedian <- rowMedians(txiTxFiltered$length[, normalIndex])
tx_normalCountMat <- rep(tx_normalCountMedian, length(singleTumorIndex))
tx_normalCountMat <- matrix(tx_normalCountMat, length(tx_normalCountMedian), length(singleTumorIndex))
colnames(tx_normalCountMat) <- pairedNormalSample
tx_normalAbundMat <- rep(tx_normalAbundMedian, length(singleTumorIndex))
tx_normalAbundMat <- matrix(tx_normalAbundMat, length(tx_normalAbundMedian), length(singleTumorIndex))
colnames(tx_normalAbundMat) <- pairedNormalSample
tx_normalLenMat <- rep(tx_normalLenMedian, length(singleTumorIndex))
tx_normalLenMat <- matrix(tx_normalLenMat, length(tx_normalLenMedian), length(singleTumorIndex))
colnames(tx_normalLenMat) <- pairedNormalSample
txiTxFiltered$abundance <- cbind(txiTxFiltered$abundance, tx_normalAbundMat)
txiTxFiltered$counts <- cbind(txiTxFiltered$counts, tx_normalCountMat)
txiTxFiltered$length <- cbind(txiTxFiltered$length, tx_normalLenMat)

ddsTx <- DESeqDataSetFromTximport(txiTxFiltered, designMat, ~age + sex + condition)
select <- rowSums(counts(ddsTx)) >=10
ddsTx <- ddsTx[select,]
ddsTx$condition <- relevel(ddsTx$condition, ref = "Normal")
ddsTx <- DESeq(ddsTx)
ddsTx <- ddsTx[which(mcols(ddsTx)$betaConv),]
#resLFC <- lfcShrink(ddsTx, coef = "condition_Tumor_vs_Normal", type = "apeglm")
resTx <- results(ddsTx)
resOrdered <- resTx[order(resTx$padj),]
write.csv(as.data.frame(resOrdered), file = "DET.00.paired_results_BY_FDR.csv")
ntdTx <- normTransform(ddsTx)
select <- which(abs(resOrdered$log2FoldChange) > 1 & resOrdered$padj < 0.01)
pheatmap(assay(ntdTx)[select, 1:58], filename = "./results/DET.Paired_cluster.pdf", show_rownames = FALSE, annotation_col = anno, fontsize_col = 6, cutree_cols = 3)

library("DEXSeq")
library("BiocParallel")
BPPARAM = SnowParam(20)
B100_CAIsoforms <- read.table(filesTx[1], header = TRUE)
dxd <- DEXSeqDataSet(countData = round(txiTxFiltered$counts), sampleData = designMat, design = ~sample + exon + condition:exon, featureID = rownames(txiTxFiltered$counts), groupID = B100_CAIsoforms$gene_id[txSelect])
levels(designMat$sample) <- sub("-",".", levels(designMat$sample))
levels(designMat$sample) <- sub("-",".", levels(designMat$sample))
rownames(designMat) <- designMat$sample
dxd <- estimateSizeFactors(dxd)
formula1 = ~sample + exon + patient:exon + condtion:exon
formula0 = ~sample + exon + patient:exon
dxd <- estimateDispersions(dxd, formula = formula1, BPPARAM = BPPARAM)
dxd <- testForDEU(dxd, fullModel = formula1, reducedModel = formula0, BPPARAM = BPPARAM)
select <- which(mcols(dxd)$fullBetaConv)
dxd <- dxd[select,]
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM = BPPARAM)
dxr <- DEXSeqResults(dxd)
qval_dtu <- perGeneQValue(dxr)
dxrSorted <- dxr[order(dxr$padj, decreasing = FALSE),]
write.csv(dxrSorted, file = "./results/DTU.00.paired_results_by_FDR.csv")
save.image("../r-image/BLCA_DTU.RData")
