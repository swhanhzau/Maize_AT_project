#use normalized expression data
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
exprMat='input.txt'
type = "unsigned"
corType = "pearson"
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T,quote="", comment="", check.names=F)
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25),na.rm=T)[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
   # Optionally, print the gene and sample names that were removed:
   if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", 
                      paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
   if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", 
                      paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
   # Remove the offending genes and samples from the data:
   dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
 }
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")

powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(dataExpr, powerVector=powers,networkType=type, verbose=5)



power = sft$powerEstimate
power
if (is.na(power)){
   power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
           ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
           ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
           ifelse(type == "unsigned", 6, 12))
           )
           )
 }
maxPOutliers = ifelse(corType=="pearson",1,0.05)
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                        TOMType = type, minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs=TRUE, corType = corType, 
                        maxPOutliers=maxPOutliers,
                        saveTOMFileBase = paste0(exprMat, ".tom"),
                        verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
   as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
cyt = exportNetworkToCytoscape(TOM,
              edgeFile = paste(exprMat, "0.15.edges.txt", sep=""),
              nodeFile = paste(exprMat, "0.15.nodes.txt", sep=""),
              weighted = TRUE, threshold = 0.15,
              nodeNames = probes, nodeAttr = moduleColors)
cyt = exportNetworkToCytoscape(TOM,
              edgeFile = paste(exprMat, "0.2.edges.txt", sep=""),
              nodeFile = paste(exprMat, "0.2.nodes.txt", sep=""),
              weighted = TRUE, threshold = 0.2,
              nodeNames = probes, nodeAttr = moduleColors)
cyt = exportNetworkToCytoscape(TOM,
              edgeFile = paste(exprMat, "0.1.edges.txt", sep=""),
              nodeFile = paste(exprMat, "0.1.nodes.txt", sep=""),
              weighted = TRUE, threshold = 0.1,
              nodeNames = probes, nodeAttr = moduleColors)
cyt = exportNetworkToCytoscape(TOM,
              edgeFile = paste(exprMat, "0.3.edges.txt", sep=""),
              nodeFile = paste(exprMat, "0.3.nodes.txt", sep=""),
              weighted = TRUE, threshold = 0.3,
              nodeNames = probes, nodeAttr = moduleColors)
cyt = exportNetworkToCytoscape(TOM,
              edgeFile = paste(exprMat, "all.edges.txt", sep=""),
              nodeFile = paste(exprMat, "all.nodes.txt", sep=""),
              weighted = TRUE, threshold = 0,
              nodeNames = probes, nodeAttr = moduleColors)
#####plot
pdf('plot1.pdf')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
pdf('plot2.pdf')
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      xlab="Soft Threshold (power)",
      ylab="Scale Free Topology Model Fit,signed R^2",type="n",
      main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
      main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
      cex=cex1, col="red")
dev.off()
pdf('plot4.pdf')
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                     "Module colors",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05)
dev.off()
pdf('plot5.pdf')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                       marDendro = c(3,3,2,4),
                       marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                       xLabelsAngle = 90)
dev.off()
pdf('plot6.pdf')
TOMplot(plotTOM, net$dendrograms, moduleColors, 
         main = "Network heatmap plot, all genes")
dev.off()
########KME
dataKME=signedKME(dataExpr, MEs, outputColumnName="kME_MM.")
write.csv(dataKME, "kME.csv")
hubs=chooseTopHubInEachModule(dataExpr,moduleColors)
HubGenes=chooseTopHubInEachModule(dataExpr,moduleColors)
write.table (HubGenes,file = "HubGenes_of_each_module.xls",quote=F,sep='\t')
connet=abs(cor(dataExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)
head(Alldegrees1)
####cytoscape
TOM = TOMsimilarityFromExpr(dataExpr, power = power)
modules = c("brown", "blue","black","green","red","turquoise","yellow")
probes = names(dataExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule])
