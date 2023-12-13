# NMR deseq
setwd('~/Projects/naked_mole_rat/')

mouse_orthology <- read.table("genome/nmr_mouse_orthologous_gene_cleaned.txt",sep="\t",header=T) 
de_result <- read.table("du/",sep=",",header=T,row.names = 1)


# change column name

de_result$nmr.gene_id <- gsub("GeneID:","",de_result$GeneName)

de_out <- merge(de_result,mouse_orthology[,-2],by="nmr.gene_id",all.x=T)

nrow(de_result)
nrow(de_out)
write.table(de_out[,-1],"du/out/all_samples_FeatureCount_count_Group3_DEgenes_P05LFC2_w_orthology.csv",quote = F,row.names = F,sep=",")
