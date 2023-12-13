# NMR deseq
setwd('~/Projects/naked_mole_rat/')

## deseq2
library(rtracklayer)
library(tidyverse)
#----------------------------------
# make count matrix
#----------------------------------
sample_list <- read.table("sample_list",sep="\t",header=F)

# get gene_id gene_name metadata
gtf_gr <- rtracklayer::import("genome/Mus_musculus.GRCm39.107.chr.gtf")
gtf_df <- as.data.frame(gtf_gr)
genes_mouse<- unique(gtf_df[,c("gene_id","gene_name")])

# get gene_id gene_name metadata
gtf_gr_refseq <- rtracklayer::import("genome/GCF_000247695.1_HetGla_female_1.0_genomic.gtf")
gtf_df_refseq <- as.data.frame(gtf_gr_refseq)
genes_refseq <- unique(gtf_df_refseq[,c("db_xref","gene")])

genes_refseq$db_xref <- gsub("GeneID:","",genes_refseq$db_xref)

colnames(genes_refseq) <- c("nmr.gene_id","nmr.symbol")
colnames(genes_mouse) <- c("mouse.gene_id","mouse.symbol")

# clean orthology table
mouse_orthology <- read.table("genome/nmr_mouse_orthologous_gene.csv",sep=",",header=T) 

mouse_orthology_1 <- merge(mouse_orthology,genes_mouse,by="mouse.gene_id")
mouse_orthology_2 <- merge(mouse_orthology_1,genes_refseq,by="nmr.gene_id")


# mouse_orthology_2 <- mouse_orthology[,c(1,6,2,5)]
# colnames(mouse_orthology_2)[c(2,4)] <- c("nmr.symbol","mouse.symbol")

write.table(mouse_orthology_2,"genome/nmr_mouse_orthologous_gene_cleaned.txt",
            quote = F,row.names = F,sep="\t")



