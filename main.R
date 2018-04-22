basePath="E:/data/tcga/lihc_tcga"
infile_path=paste(basePath,"lihc_tcga",sep="/")
shell.exec(infile_path)

##读取数据ran_seq_median
infile=list.files(path = infile_path,pattern ="data_RNA_Seq_v2_expression_median" ,full.names = TRUE)
RNAseq=read.csv(file = infile,sep = "\t")
fix(RNAseq)
RNAseq1=RNAseq[,-1]
rownames(RNAseq1)=RNAseq1$Entrez_Gene_Id
RNAseq1=RNAseq1[,-1]
fix(RNAseq1)
##结果输出的位置
out_put_dir=paste(basePath,"WGCNA",sep="/")
shell.exec(out_put_dir)

library(parallel)
# 
# library(limma)
# RNAseq_voom = voom(RNAseq1)$E
#transpose matrix to correlate genes in the following
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq1,1,mad), decreasing = T)[1:5000],])


#similarity measure between gene profiles: biweight midcorrelation
library(WGCNA)
s = abs(bicor(WGCNA_matrix))


library(WGCNA)
WGCNAnThreads(2)

