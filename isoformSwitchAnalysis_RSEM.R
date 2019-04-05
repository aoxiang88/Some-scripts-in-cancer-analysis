library("IsoformSwitchAnalyzeR")
library("BSgenome.Hsapiens.NCBI.GRCh38")
setwd("/home/aoxiang/data/bladder-cancer/rsem_results")
sampleID <- read.table("../clean_fq/bam_58.list")
tid <- read.table("./isoform_switch_prefilt/transcript_id.list", header=TRUE)
files <- file.path(getwd(), paste0("/expression_data/", sampleID$V1,".isoforms.results"))
tx <- tximport::tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)
deGeneID <- read.table("DESeq/DESeq_DE_geneID_sig_ensemblID.txt", header = TRUE)
pairGeneID <- read.table("./expression_data/gene_id_ensembl.txt", header = TRUE)
select <- which(! pairGeneID$ensemblID %in% deGeneID$geneID )
tx_count <- tx$counts[select, ]
tx_count <- as.data.frame(tx_count)
rownames(tx_count) <- tid$transcript_id[select]
colnames(tx_count) <- sampleID$V1
tx_count$isoform_id <- rownames(tx_count)
clinical <- read.table("../clean_fq/clinical_characteristics_RNAseq.list", header=TRUE)
design_mat <- data.frame(sampleID = sampleID$V1, condition = clinical$State, sex = clinical$Sex, age = clinical$Age, sample = clinical$Sample)
switch_list <- importRdata(tx_count, designMatrix = design_mat, isoformExonAnnoation = "/home/aoxiang/library/Homo_Sapiens.GRCh38.91.primary.gtf")
switch_list_filtered <- preFilter(switch_list, geneExpressionCutoff = 0.5, isoformExpressionCutoff = 0, removeSingleIsoformGenes = TRUE)
switch_list_analyzed <- isoformSwitchTest(switch_list_filtered)
switch_list_analyzed <- analyzeORF(switch_list_analyzed, genomeObject = Hsapiens)
switch_list_analyzed <- extractSequence(switch_list_analyzed, genomeObject = Hsapiens, pathToOutput = './fasta')
# signalP(wd:~/data/bladder-cancer/rsem_results/fasta): $signalp -f summary isoformSwitchAnalyzeR_isoform_AA.fasta > signalp_result.out &
# cpat(wd:~/data/bladder-cancer/rsem_results/fasta): $cpat.py -g isoformSwitchAnalyzeR_isoform_nt.fasta -d Human_logitModel.RData -x Human_Hexamer.tsv -o cpat_result.out&
# Pfam(wd: ~/data/bladder-cancer/rsem_results/fasta): $pfam_scan.pl -fasta isoformSwitchAnalyzeR_isoform_AA.fasta -dir ~/_bio_Tools/PfamScan -outfile pfam_result.out -cpu 20&
switch_list_analyzed <- analyzePFAM(switchAnalyzeRlist = switch_list_analyzed, pathToPFAMresultFile = "./fasta/pfam_result.out", showProgress = TRUE)
switch_list_analyzed <- analyzeSignalP(switchAnalyzeRlist = switch_list_analyzed, pathToSignalPresultFile = "./fasta/signalp_result.out")
switch_list_analyzed <- analyzeCPAT(switchAnalyzeRlist = switch_list_analyzed, pathToCPATresultFile = "./fasta/cpat_result.out", codingCutoff = 0.725, removeNoncodinORFs = TRUE)
switch_list_analyzed <- analyzeIntronRetention(switch_list_analyzed)

consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')
switch_list_analyzed <- analyzeSwitchConsequences(switch_list_analyzed, consequencesToAnalyze = consequencesOfInterest, dIFcutoff = 0.2)
