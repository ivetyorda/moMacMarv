
#set directory to location of TSV values----
setwd("C:/Users/user/folder")
files = list.files()

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

#reorder files
files = files[c(1,2,3,4)]
list.files()
all(file.exists(files))

#use these to set sample names in order 
samples=c("sample1", "sample2", "sample3", "sample4")

# Load the GTF/GFF file gtf_file
library(rtracklayer)
gtf_file <- "C:/Users/user/folder/GCF_000857325.2_ViralProj15199_genomic.gtf"
gtf <- import(gtf_file) 

# Extract transcript-to-gene mappings
tx2gene <- gtf[gtf$type == "transcript"] 
tx2gene <- data.frame(transcript_id = tx2gene$transcript_id, gene_id = tx2gene$gene_id)

#import transcript abundance data # removed "ignoreTxVersion = TRUE" which led to an error in annotation, in this case, there are multiple versions for each transcript.
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

#set up metadata----
cell_type = as.factor(c(rep("treatment1",2), rep("treatment2",2)))
cell_type = data.frame(cell_type)

#add metadata columns where a given cell type has a 1 and all other cell types have a 0

#Here, two cell types (or comparison groups) are defined, mac_control and dend_control. 
comparisons = c("treatment1", "treatment2") 

#A data frame is created, meta_data, with two columns: sample and cell_type. It is assumed that samples and cell_type are already defined elsewhere in your code.
meta_data = data.frame(sample = samples, cell_type = cell_type)

meta_data[comparisons] = 0
for (i in 1:length(comparisons)){
  meta_data[which(meta_data$cell_type==comparisons[i]),comparisons[i]]=1
}

meta_data[comparisons]=lapply(meta_data[comparisons], factor)

#general dds----
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta_data,
                                design = ~ cell_type)
as.data.frame(colData(dds))

dds<- DESeq(dds)
class(dds)

rld<-rlog(dds)

#save normalized gene counts
write.csv(as.data.frame(assay(rld)), file="filename.csv" )

