# NMR deseq
setwd('~/Projects/naked_mole_rat/human_subject')

## deseq2
library(tximport)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggvenn)
library(affy)
library(limma)
library(hgu133plus2.db)


phenoData <- read.AnnotatedDataFrame("sample_meta_ALKvsNormal.txt")
celpath <- "CEL_ALK_normal"
cel_file <- ReadAffy(celfile.path = celpath, phenoData = phenoData)
eset <- rma(cel_file,phenoData=phenoData)

# annotation
exp_df <- data.frame(exprs(eset))
annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "),
                    SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "),
                    DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))

all <- merge(annot,exp_df,by.x=0,by.y=0,all=T)

combn <- factor(pData(phenoData)[,1],levels=c("ALK","Normal"))
design <- model.matrix(~combn)
colnames(design) <- c("ALK","ALKvsNormal")

fit <- lmFit(eset, design)

efit <- eBayes(fit)
topTable(efit,coef="ALKvsNormal",adjust="BH") 

out_df <- topTable(efit,coef="ALKvsNormal",adjust="BH",n=Inf) 

# add gene name

out_final <- merge(out_df,annot,by.x=0,by.y=0,all=T)

write.table(out_final,"human_exp_results_ALK_normal.txt",sep="\t",quote = F,row.names = F)


