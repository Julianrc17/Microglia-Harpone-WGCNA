joint<-read.table("combatExpr.txt", header =T)
head(joint)

samplesdf<-read.table('joint_samplesdf.txt', header =T, sep="\t")
head(samplesdf)


samplesdf_filtered <- subset(samplesdf, Genotype %in% c("E2WT", "E3WT", "E4WT", "KO"))
str(samplesdf)
str(samplesdf_filtered)
annotation<-read.table("GO_screen_annot.txt", header=T, sep="\t")
head(annotation)

microsamples<-samplesdf_filtered[which(samplesdf_filtered$Cell_type=="microglia"),]
microdf<-joint[,which(colnames(joint) %in% rownames(microsamples))]

colnames(microdf)
rownames(microsamples)

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = microdf;
# Take a quick look at what is in the data set:
dim(femData);
names(femData);

femData0<-merge(femData,annotation,by.x="row.names",by.y="transcript_cluster_id")
head(femData0)
rownames(femData0)<-make.names(femData0$SYMBOL,unique=T)
femData0 <- femData0[, !(colnames(femData0) %in% c("gene", "SYMBOL", "Row.names"))]
head(femData0)


datExpr0 <- as.data.frame(t(femData0), stringsAsFactors = FALSE)

datExpr0[] <- lapply(datExpr0, function(x) as.numeric(as.character(x)))

datExpr0 <- datExpr0[rowSums(is.na(datExpr0)) < ncol(datExpr0), ]  # Eliminar muestras con solo NA
datExpr0 <- datExpr0[, colSums(is.na(datExpr0)) < nrow(datExpr0)]  # Eliminar genes con solo NA

str(datExpr0)
head(datExpr0)
                     
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2.5,
cex.axis = 1.5, cex.main = 2)



# Plot a line to show the cut
# abline(h = 80, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

traitData = samplesdf_filtered
traitData$KO <- ifelse(traitData$Genotype == "KO", 1, 0)
traitData$E2WT <- ifelse(traitData$Genotype == "E2WT", 1, 0)
traitData$E3WT <- ifelse(traitData$Genotype == "E3WT", 1, 0)
traitData$E4WT <- ifelse(traitData$Genotype == "E4WT", 1, 0)

datTraits = traitData[rownames(traitData) %in% rownames(datExpr), ]

datTraits <- datTraits[, c("KO", "E2WT", "E3WT", "E4WT")]

print(head(datTraits))
unique(traitData$Genotype)
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Microglia-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr, power = 12,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "femaleMouseTOM",
verbose = 3)

mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "Microglia-02-networkConstruction-auto.RData")

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Microglia-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Microglia-02-networkConstruction-auto.RData");
lnames
