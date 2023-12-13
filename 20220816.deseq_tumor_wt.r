# NMR deseq
setwd('~/Projects/naked_mole_rat/')

## deseq2
library(tximport)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggvenn)

raw_count <- read.table("du/all_samples_FeatureCount_count_read_raw_counts_geneInfo.csv",sep=",",header=T)
raw_count <- raw_count[which(raw_count$GeneType=="protein_coding"),]
id_gene <- raw_count[,c(2,4)]

row.names(raw_count) <- raw_count[,2]
raw_count <- raw_count[,-c(1:4)]
sample_data <- data.frame(row.names = colnames(raw_count),
                          genotype = factor(c(rep("Tumor",11),"Control"),
                                            levels=c("Control","Tumor")))
sample_data$name <- paste(sample_data$genotype,
                          c("1","2","3","4","5","6","7","8","9","10","11",
                            "1"),sep="_")

#===============================================================
# 1. gene expression changes between the Tumor and Control cells
#===============================================================

dds <- DESeqDataSetFromMatrix(countData = raw_count,
                              colData = sample_data,
                              design = ~ genotype)
dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds)

# lfc shrinking
lfcs_res<- lfcShrink(dds,coef="genotype_Tumor_vs_Control", res=res)

# add mouse gene
mouse_orthology <- read.table("genome/nmr_mouse_orthologous_gene_cleaned.txt",sep="\t",header=T) 

colnames(mouse_orthology)[c(2,4)] <-c("Symbol","mouse_gene_name")
mouse_orthology$Symbol <- toupper(mouse_orthology$Symbol)
mouse_orthology$mouse_gene_name <- toupper(mouse_orthology$mouse_gene_name)

colnames(id_gene) <- c("Symbol","nmr.gene_id")
  
id_gene$nmr.gene_id <- gsub("GeneID:","",id_gene$nmr.gene_id)

lfcs_res$Symbol <- row.names(lfcs_res)
out_df <- merge(as.data.frame(lfcs_res),id_gene,by="Symbol")

out_df_2 <- merge(out_df,mouse_orthology,by="nmr.gene_id")

write.table(out_df_2,"results/tumor_control_DEgene.txt",sep="\t",quote = F,row.names = F)
#================================================
# pca & correlation
#================================================
library(corpcor)

# vsd_nf2null <- counts(dds, normalized=TRUE)
# vsd_nf2null <- assay(rlog(dds, blind=FALSE))
vsd <- assay(rlog(dds, blind=FALSE))
svd_tbl <- wt.scale(t(vsd),center=TRUE,scale=TRUE)
svd_run <- fast.svd(svd_tbl)
ds <- round(svd_run$d^2/sum(svd_run$d^2)*100,2)

# PCA
pdf("results/PCA_tumor_control.pdf", width = 8,height = 6)
par(mfrow=c(1,2),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
plot(svd_run$u[,1],svd_run$u[,2],col=as.integer(as.factor(sample_data[,1]))+1,pch=20,
     xlab=paste0('PC1: ',ds[1],'%'),ylab=paste0('PC2: ',ds[2],'%'),main='')
text(svd_run$u[,1],svd_run$u[,2],sample_data[,2],cex=0.6,pos=1)
#legend('topright',levels(as.factor(sample_data[,1])),col=as.integer(as.factor(levels(as.factor(sample_data[,1]))))+1,pch=20,bty='n',cex=0.6)
plot(svd_run$u[,1],svd_run$u[,3],col=as.integer(as.factor(sample_data[,1]))+1,pch=20,
     xlab=paste0('PC1: ',ds[1],'%'),ylab=paste0('PC3: ',ds[3],'%'),main='')
text(svd_run$u[,1],svd_run$u[,3],sample_data[,2],cex=0.6,pos=1)
dev.off()

# correlation
pdf("results/correlation_tumor_control.pdf", width = 8,height = 4)
par(mfrow=c(1,2),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
for(i in 1:2){
    for(j in 1:1){
        x <- (i-1)*2+j
        y <- (i-1)*2+j+1
        if(y > i*2) y <- y-2
        cat(x,y,'\n')
        #plot(vsd[,x],vsd[,y],pch=20,col=rgb(0,0,0,alpha=0.1),xlab=sample_data$name[x],ylab=sample_data$name[y])
        plot(log2(assay(dds)[,x]+1),log2(assay(dds)[,y]+1),pch=20,col=rgb(0,0,0,alpha=0.1),
             xlab=sample_data$name[x],
             ylab=sample_data$name[y])
        abline(a=0,b=1)
    }
}
dev.off()

pdf("results/hist_tumor_control.pdf", width = 8,height = 8)
par(mfrow=c(4,3),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
for(i in 1:12){
    hist(vsd[,i],n=100,xlab=sample_data$name[i],main='')
}
dev.off()

#================================================
# make result
#================================================
# lfcs list
lfcs_list <- list(#NF2null = lfcs_NF2null,
                  WT=lfcs_res)
                  #interaction=lfcs_interaction,
                  #NF2null_WT_control=lfcs_NF2null_WT_control,
                  #NF2null_WT_JQ1=lfcs_NF2null_WT_JQ1)

# DEseq2 normalized count
dat_normal<- counts(dds, normalized=TRUE)
dat <- cbind(lfcs_res,dat_normal)

pdf("results/MA_tumor_controal.pdf", width = 4,height = 4)
par(mfrow=c(1,1),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
plotMA(lfcs_res)

# for(i in 1:length(lfcs_list)){
#     DESeq2::plotMA(lfcs_list[[i]], ylim=c(-5,5),main=names(lfcs_list)[i])
#     tmp <- data.frame(lfcs_list[[i]][,c(1,2,4,5)],stringsAsFactors=F)
#     colnames(tmp) <- paste0(names(lfcs_list)[i],c('.mean','.lfc','.pv','.fdr'))
#     #dat <- tmp[match(rownames(dat),rownames(tmp)),]
#     dat <- cbind(dat,tmp[match(rownames(dat),rownames(tmp)),])
# }
dev.off()

# add gene names
dat$Symbol <- rownames(dat)
dat <- merge(id_gene,dat,by="gene_id")

write.table(dat,file='results/result_normalized_tumor_control.txt',sep='\t',row.names=F,quote=F)

# idx <- genes$gene_type[match(rownames(dat),genes$gene_id)]=='protein_coding' ## only protein coding genes
# plot(log2(dat[idx,4]+1),log2(dat[idx,5]+1),pch=20,col=rgb(0,0,0,alpha=0.1))

#------------------------------------------------------------------------------
# volcano plots
#------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)

dat <- read.table("results/result_normalized_tumor_control.txt",sep='\t',header=T)
mouse_orthology <- read.table("genome/nmr_mouse_orthologous_gene_cleaned.txt",sep="\t",header=T) 

colnames(mouse_orthology)[c(2,4)] <-c("Symbol","mouse_gene_name")
mouse_orthology$Symbol <- toupper(mouse_orthology$Symbol)
mouse_orthology$mouse_gene_name <- toupper(mouse_orthology$mouse_gene_name)

dat_orthology_tmp <- merge(dat,mouse_orthology[,-1],by="gene_name",all=T)
dat_orthology <- dat_orthology_tmp[which(dat_orthology_tmp$gene_id %in% dat$gene_id),]

#write.table(dat_orthology,"results/result_normalized_with_mouse_orthology_gene.txt",row.names = F,quote = F,sep="\t")

# one condition
dat_chosen <- dat_orthology[!is.na(dat_orthology$mouse_gene_name),]
vol_df <- dat_chosen[,c(13,3,4,7)]
colnames(vol_df) <- c("Symbol","mean","lfc","fdr")

# remove NA
vol_df <- vol_df[which(!is.na(vol_df$fdr) | !is.na(vol_df$fdr)),]
vol_df$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
vol_df$diffexpressed[vol_df$lfc > 0.8 & vol_df$fdr < 0.05 ] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
vol_df$diffexpressed[vol_df$lfc < -0.8 & vol_df$fdr < 0.05] <- "DOWN"

# get significant (lfc > 0.8, fdr < 0.05, mean > 50)
vol_df_sig <- vol_df[which(vol_df$fdr < 0.05 & abs(vol_df$lfc) > 0.58 & vol_df$mean > 50),]

# get top 50
vol_df_top <- head(vol_df_sig[order(-abs(vol_df_sig$lfc)),],50)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


volcano <- ggplot(data=vol_df,aes(x=lfc,y=-log10(fdr),
                       col = diffexpressed, label=Symbol))+
    geom_point(size=1)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size = 12))+
    geom_text_repel(data = vol_df_top,
                    max.overlaps = 50,size = 3)+
    geom_vline(xintercept=c(-0.8, 0.8), col="red",linetype='dashed') +
    geom_hline(yintercept=-log10(0.05), col="red",linetype='dashed')+
    scale_color_manual(values = mycolors)+
    ggtitle("WT vs ALK+")+
    ylab("-log10(FDR)")+
    xlab("log2 fold change")

pdf(file = "results/volcano_2.pdf",width = 8, height = 7)
volcano
dev.off()


#------------------------------------------------------------------------------
# heatmap
#------------------------------------------------------------------------------
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendextend)

# one condition
temp_dat <- dat_chosen[,c(13,8:11,3,4,7)]
colnames(temp_dat)[2:8] <- c("cat_372.1_S5","cat_372.2_S6","cat_372.3_S7",
                              "cat_372.4_S8","mean","lfc","fdr")
# temp_dat <- dat[,c(1:2,
#                    6:8,
#                    12:14,
#                    33,34,36)]
#colnames(temp_dat)[9:11] <- c("mean","lfc","fdr")
sig_gene <- temp_dat[which(temp_dat$fdr < 0.05 &
                               abs(temp_dat$lfc)>0.58 &
                               temp_dat$mean> 50),]

write.table(sig_gene,"results/sig_diff_gene.txt",quote = F,row.names = F,sep = "\t")

# get top 50
sig_gene_top <- head(sig_gene[order(-abs(sig_gene$lfc)),],50)

# make heatmap matrix
sig_gene_plot_df <- sig_gene_top[,c(1,2:5)]
#sig_gene_plot_df <- sig_gene_top[,c(2:14)]
row.names(sig_gene_plot_df) <- sig_gene_plot_df[,1]

# make sample groups
sample_group <- data.frame(Conditions=rep(c("WT",
                                            "ALK"),
                                          c(2,2)))
row.names(sample_group) <- colnames(sig_gene_plot_df[,-1])
sample_group$Conditions <- factor(sample_group$Conditions,levels = c("WT",
                                                                     "ALK"))
group_color <- brewer.pal(4,"Dark2")
annotation_color <- list(Conditions = c(WT=group_color[1],
                                        ALK=group_color[2]))

heatmap<- pheatmap(sig_gene_plot_df[,-1],
         #color = rev(brewer.pal(9,"RdBu")),
         border_color = "black",
         show_rownames = T,
         #cluster_cols = F,
         #cellwidth = 15,
         #cellheight = 1,
         labels_row = sig_gene_plot_df$gene_id,
         annotation_colors = annotation_color,
         #gaps_col = 3,
         annotation_col = sample_group,
         scale="row",
         cellwidth=30,
         main = "WT vs. ALK+ (Top 50)",
         treeheight_col = 0
         #fontsize_row = 0.4,
         #annotation_row = gene_col
         )

pdf(file = "results/heatmap_2.pdf",width = 8, height = 8)
heatmap
dev.off()

# vinn diagram
main_result <- read.table("results_v2/main_result_normalized.txt",sep="\t",header=T)

control_top_50 <- read.table("results_v2/NF2null_WT_control_top50",sep="\t",header=F)
jq_top_50 <- read.table("results_v2/NF2null_WT_JQ1_top50",sep="\t",header=F)

p1_plot_df <- list(Control = control_top_50$V1,
                   JQ1 = jq_top_50$V1)

p1 <- ggvenn(p1_plot_df)+
    ggtitle("NF2null vs. WT")

pdf(file = "results_v2/JQ1_control_vinn.pdf",width = 4, height = 4)
p1
dev.off()
