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
WGCNA_matrix = t(RNAseq1[order(apply(RNAseq1,1,mad), decreasing = T)[1:10000],])


#similarity measure between gene profiles: biweight midcorrelation
library(WGCNA)
s = abs(bicor(WGCNA_matrix))


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

#calculation of adjacency matrix
beta = 3
a = s^beta
rm(a);rm(s);
#dissimilarity measure
w = 1-a

#create gene tree by average linkage hierarchical clustering 
geneTree = hclust(as.dist(w), method = 'average');rm(w)
#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
#assign module colours
module.colours = labels2colors(modules)

#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')



library(ape)
#calculate eigengenes
MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours, excludeGrey = FALSE)$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');

#plot the result with phytools package
par(mar=c(2,2,2,2))
a=plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), bg = levels(as.factor(module.colours)))

