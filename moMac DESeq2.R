#set directory to location of TSV values----
setwd("C:/Users/user/folder")
files = list.files()

#reorder files
files = files[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)]
list.files()
all(file.exists(files))

#use these to set sample names in order 
samples=c("control1", "control2", "control3", "control4", "control5", "LPS1", "LPS2", "LPS3", "LPS4", "LPS5", "MARV1", "MARV2", "MARV3", "MARV4", "MARV5")

#open packages
library(tximport)
library(readr)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
library(apeglm)
library("BiocParallel")
register(SnowParam(4))
library("pheatmap")
library("tidyverse")
library(edgeR)
library(RColorBrewer)

# Load the GTF/GFF file gtf_file
library(rtracklayer)
gtf_file <- "C:/Users/user/folder/GCF_000001405.40_GRCh38.p14_genomic.gtf"
gtf <- import(gtf_file) 

# Extract transcript-to-gene mappings

tx2gene <- gtf[gtf$type == "transcript"] 
tx2gene <- data.frame(transcript_id = tx2gene$transcript_id, gene_id = tx2gene$gene_id)

#import transcript abundance data 
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

#set up metadata----
cell_type = as.factor(c(rep("control",5), rep("lps",5), rep("marv",5)))
cell_type = data.frame(cell_type)

#add metadata columns where a given cell type has a 1 and all other cell types have a 0

#Here, two cell types (or comparison groups) are defined, mac_control and dend_control. 
comparisons = c("control", "lps")

#A data frame is created, meta_data, with two columns: sample and cell_type. It is assumed that samples and cell_type are already defined elsewhere in your code.
meta_data = data.frame(sample = samples, cell_type = cell_type)

meta_data[comparisons] = 0 
for (i in 1:length(comparisons)){
  meta_data[which(meta_data$cell_type==comparisons[i]),comparisons[i]]=1
}

meta_data[comparisons]=lapply(meta_data[comparisons], factor)
meta_data

#general dds----
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta_data,
                                design = ~ cell_type)
as.data.frame(colData(dds))

dds<- DESeq(dds)
class(dds)

rld<-rlog(dds)
write.csv(as.data.frame(assay(rld)), file="filename.csv" )

#make DESeqDataSet for control comparison----
dds_control = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ control)

#run differential expression analysis with default parameters
dds_control = DESeq(dds_control)
resultsNames(dds_control)
res_control = results(dds_control, contrast = c("control", "1", "0"), name="control_1_vs_0")
?results
res_control
write.csv( as.data.frame(res_control), "filename_control.csv" )

#make DESeqDataSet for lps comparison----
dds_lps = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ lps)

#run differential expression analysis with default parameters
dds_lps = DESeq(dds_lps)
resultsNames(dds_lps)
res_lps = results(dds_lps, contrast = c("lps", "1", "0"), name="lps_1_vs_0")
res_lps
write.csv( as.data.frame(res_lps), "filename_lps.csv" )

#filter and save only significantly-differentially-expressed genes
res_lps_Sig <- subset(res_lps, padj < 0.05)
res_lps_Sig
write.csv( as.data.frame(res_lps_Sig), "filename_lps_sig.csv" )

#Start second comparison
comparisons = c("control", "marv")

meta_data2 = data.frame(sample = samples, cell_type = cell_type)

meta_data2[comparisons] = 0
for (i in 1:length(comparisons)){
  meta_data2[which(meta_data2$cell_type==comparisons[i]),comparisons[i]]=1
}


meta_data2[comparisons]=lapply(meta_data2[comparisons], factor)
meta_data2

#make DESeqDataSet for marv comparison----
dds_marv = DESeqDataSetFromTximport(txi, colData = meta_data2, design = ~ marv)

#run differential expression analysis with default parameters
dds_marv = DESeq(dds_marv)
resultsNames(dds_marv)
res_marv = results(dds_marv, contrast = c("marv", "1", "0"), name="marv_1_vs_0")
res_marv
write.csv( as.data.frame(res_marv), "filename_marv.csv" )
res_marv_Sig <- subset(res_marv, padj < 0.05)
res_marv_Sig
write.csv( as.data.frame(res_marv_Sig), "filename_marv_sig.csv" )

#make volcano plot for LPS vs mock. Then repeat code for MARV vs mock

library(ggplot2)
library(plotly)

topT <- as.data.frame(res_lps)

#rename column 1 name

topT2 <- cbind(rownames(topT), topT)
rownames(topT2) <- NULL
colnames(topT2) <- c("geneID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")

highlight_df <-topT2 %>% dplyr::filter(log2FoldChange>=2)
highlight_df_low <-topT2 %>% dplyr::filter(log2FoldChange<=-2)
highlight_df_mid <-topT2 %>% dplyr::filter(log2FoldChange>=-2 & log2FoldChange<=2)

any(is.na(topT2))

#If the answer was TRUE, then you have to delete NA's as below:

topT2 = topT2[-which(is.na(topT2))]

vplot <- ggplot(topT2) +
  aes(y=-log10(padj), x=log2FoldChange, text = paste("Symbol:", geneID)) +
  geom_point(data=highlight_df_mid, aes(x=log2FoldChange,y=-log10(padj)),alpha=0.5, color="#CCCCCC",size=4)+ scale_y_continuous(limits=c(0,300)) + scale_x_continuous(limits=c(-12,13)) +
  geom_point(data=highlight_df, aes(x=log2FoldChange,y=-log10(padj)), alpha=0.5,color="turquoise4",size=4) +
  geom_point(data=highlight_df_low, aes(x=log2FoldChange,y=-log10(padj)),alpha=0.5, color='firebrick',size=4) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="black", linewidth=0.8) +
  geom_vline(xintercept = 3, linetype="longdash", colour="black", linewidth=0.8) +
  geom_vline(xintercept = -3, linetype="longdash", colour="black", linewidth=0.8) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="responders vs non-responders",
       subtitle ="",
       caption=paste0("produced on ", Sys.time())) +
  theme_classic() +
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 15)+
          theme(data.title = element_text(size = 15)))
vplot
ggplotly(vplot)
ggsave("/Users/user/folder/filename.tiff", width=5, height=5)

#make a global heatmap

rowData(dds)

#filter the data to remove missing genes and non-significant genes based on Walt statistical tests

ddssubset <- subset(dds, baseMean>0)
ddssubset2 <- subset(ddssubset, WaldPvalue_cell_type_mock_vs_lps<0.05)
ddssubset3 <- subset(ddssubset, WaldPvalue_cell_type_marv_vs_lps<0.05)

rowData(ddssubset3)

rownames(ddssubset2) <- TRUE
colnames(ddssubset2) <- TRUE
pheatmap(assay(ddssubset3), scale="row", cluster_cols=FALSE, color=colorRampPalette(c("turquoise4", "white", "firebrick"))(50))





