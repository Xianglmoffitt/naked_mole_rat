# plot scripts

setwd("~/Projects/naked_mole_rat/results/GSEA/")

library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(scales)
library(data.table)

gsea_pos <- read.table("NMR_tumor_control/gsea_report_for_na_pos_1680803923881.tsv",header=T,sep = "\t")
gsea_neg <- read.table("NMR_tumor_control/gsea_report_for_na_neg_1680803923881.tsv",header=T,sep = "\t")

pos_temp <- gsea_pos[,c(1,8)]
pos_temp$score <- -log10(pos_temp$FDR.q.val)
pos_temp$group <- "Up"

neg_temp <- gsea_neg[,c(1,8)]
neg_temp$score <- log10(neg_temp$FDR.q.val)
neg_temp$group <- "Down"

barplot_df <- rbind(pos_temp,neg_temp)

barplot_df$score[which(barplot_df$score=="-Inf")] <- -3

bar_plot <- ggplot(barplot_df,aes(x = reorder(NAME, score), y = score,fill=group)) +
  # set overall appearance of the plot
  theme_classic() +
  geom_bar(stat="identity") +
  # Set main and axis titles
  ggtitle("NMR (Tumor vs. Control)") +
  xlab("GSEA mouse Hallmark gene set") +
  ylab("Negative Log10 FDR") +
  # Add a line showing the alpha = 0.05 level
  geom_hline(yintercept = -log10(0.05), size = 1, color = "black",linetype="dashed") +
  geom_hline(yintercept = log10(0.05), size = 1, color = "black",linetype="dashed") +
  # Flip the x and y axes
  coord_flip()+
  scale_fill_manual(values=c("blue","red"))+
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10))

pdf(file = "NMR_gsea_Tumor_control.pdf",width = 8, height = 8)
bar_plot
dev.off()

