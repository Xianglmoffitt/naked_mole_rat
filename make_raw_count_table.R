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
gtf_gr <- rtracklayer::import("genome/Heterocephalus_glaber_female.HetGla_female_1.0.107.gtf")
gtf_df <- as.data.frame(gtf_gr)
genes <- unique(gtf_df[,c("gene_id","gene_name")])

# replace NA gene_name to gene_id
na_idx <- which(is.na(genes$gene_name))
genes$gene_name[na_idx] <- genes$gene_id[na_idx]


# make count matrix
count_list <- list()
for (i in 1:nrow(sample_list)) {
  df_tmp <- read.table(paste0("feature_count/",sample_list$V1[i]),sep="\t",header=T)
  colnames(df_tmp)[c(1,7)] <- c("gene_id",sample_list$V1[i]) 
  count_list[[i]] <- df_tmp[,c(1,7)]
}

raw_count <- count_list %>% reduce(inner_join, by = "gene_id")

# add gene symbol
raw_count <- merge(raw_count,genes,by="gene_id")

raw_count_final <- raw_count[,c(1,6,2:5)]

write.table(raw_count_final,"results/raw_count.txt",quote=F,row.names = F,sep="\t")


