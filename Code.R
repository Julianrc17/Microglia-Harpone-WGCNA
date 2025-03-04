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

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.8,
zlim = c(-1,1),
main = paste("Module-trait relationships"))


head(moduleTraitCor)
head(moduleTraitPvalue)
# head(textMatrix)

head(cbind(data.frame(moduleTraitCor[,1]),data.frame(moduleTraitPvalue[,1])))

dfcor<-cbind(data.frame(moduleTraitCor[,1]),data.frame(moduleTraitPvalue[,1]))
colnames(dfcor)<-c("R_KO","P_KO")
head(dfcor)
write.table(dfcor,"module.correlations.txt",sep="\t")

dfcor[order(dfcor$P_KO,decreasing=F),]

# Define variable KO containing the KO column of datTrait
KO = as.data.frame(datTraits$KO);
names(KO) = "KO"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, KO, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(KO), sep="");
names(GSPvalue) = paste("p.GS.", names(KO), sep="");

# Definir los genotipos a partir de la variable datTraits
KO = as.data.frame(datTraits$KO)
E2WT = as.data.frame(datTraits$E2WT)
E3WT = as.data.frame(datTraits$E3WT)
E4WT = as.data.frame(datTraits$E4WT)

# Calcular la correlación de los genotipos con los módulos
# Para KO
geneTraitSignificance_KO = as.data.frame(cor(datExpr, KO, use = "p"))
GSPvalue_KO = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_KO), nSamples))
names(geneTraitSignificance_KO) = paste("GS.", names(KO), sep="")
names(GSPvalue_KO) = paste("p.GS.", names(KO), sep="")

# Para E2WT
geneTraitSignificance_E2WT = as.data.frame(cor(datExpr, E2WT, use = "p"))
GSPvalue_E2WT = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_E2WT), nSamples))
names(geneTraitSignificance_E2WT) = paste("GS.", names(E2WT), sep="")
names(GSPvalue_E2WT) = paste("p.GS.", names(E2WT), sep="")

# Para E3WT
geneTraitSignificance_E3WT = as.data.frame(cor(datExpr, E3WT, use = "p"))
GSPvalue_E3WT = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_E3WT), nSamples))
names(geneTraitSignificance_E3WT) = paste("GS.", names(E3WT), sep="")
names(GSPvalue_E3WT) = paste("p.GS.", names(E3WT), sep="")

# Para E4WT
geneTraitSignificance_E4WT = as.data.frame(cor(datExpr, E4WT, use = "p"))
GSPvalue_E4WT = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_E4WT), nSamples))
names(geneTraitSignificance_E4WT) = paste("GS.", names(E4WT), sep="")
names(GSPvalue_E4WT) = paste("p.GS.", names(E4WT), sep="")

# Definir los nombres de los módulos (modNames) y los colores de los módulos (moduleColors)
modNames = substring(names(MEs), 3)
moduleColors = as.vector(moduleColors)

# Graficar las correlaciones para cada módulo y cada genotipo
modules = c("red", "purple", "black", "tan", "blue", "green", "brown", "magenta")

# Función para hacer el gráfico de dispersión para cada genotipo y módulo
plot_gene_module_correlation <- function(module, genotipo, geneModuleMembership, geneTraitSignificance, moduleColors) {
  column = match(module, modNames)
  moduleGenes = moduleColors == module
  
  par(mfrow = c(1, 1))
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = paste("Gene significance for", genotipo),
    main = paste("Module membership vs. gene significance for", genotipo),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
  )
}

# Correlaciones y gráficos para cada genotipo y módulo
for (module in modules) {
  # Para KO
  plot_gene_module_correlation(module, "KO", geneModuleMembership, geneTraitSignificance_KO, moduleColors)
  
  # Para E2WT
  plot_gene_module_correlation(module, "E2WT", geneModuleMembership, geneTraitSignificance_E2WT, moduleColors)
  
  # Para E3WT
  plot_gene_module_correlation(module, "E3WT", geneModuleMembership, geneTraitSignificance_E3WT, moduleColors)
  
  # Para E4WT
  plot_gene_module_correlation(module, "E4WT", geneModuleMembership, geneTraitSignificance_E4WT, moduleColors)
}

# Exportar los resultados de las correlaciones y p-valores

# Para KO
dfcor_KO <- cbind(data.frame(geneTraitSignificance_KO[, 1]), data.frame(GSPvalue_KO[, 1]))
colnames(dfcor_KO) <- c("R_KO", "P_KO")
write.table(dfcor_KO, "module_correlations_KO.txt", sep = "\t")

# Para E2WT
dfcor_E2WT <- cbind(data.frame(geneTraitSignificance_E2WT[, 1]), data.frame(GSPvalue_E2WT[, 1]))
colnames(dfcor_E2WT) <- c("R_E2WT", "P_E2WT")
write.table(dfcor_E2WT, "module_correlations_E2WT.txt", sep = "\t")

# Para E3WT
dfcor_E3WT <- cbind(data.frame(geneTraitSignificance_E3WT[, 1]), data.frame(GSPvalue_E3WT[, 1]))
colnames(dfcor_E3WT) <- c("R_E3WT", "P_E3WT")
write.table(dfcor_E3WT, "module_correlations_E3WT.txt", sep = "\t")

# Para E4WT
dfcor_E4WT <- cbind(data.frame(geneTraitSignificance_E4WT[, 1]), data.frame(GSPvalue_E4WT[, 1]))
colnames(dfcor_E4WT) <- c("R_E4WT", "P_E4WT")
write.table(dfcor_E4WT, "module_correlations_E4WT.txt", sep = "\t")

# Crear el data frame inicial
geneInfo0 = data.frame(
  moduleColor = moduleColors,
  geneTraitSignificance,  # Ya contiene los valores para KO
  GSPvalue                # Ya contiene los valores p para KO
)

# Añadir las significancias de los otros fenotipos (E2WT, E3WT, E4WT)
geneTraitSignificance_E2WT = as.data.frame(cor(datExpr, datTraits$E2WT, use = "p"))
GSPvalue_E2WT = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_E2WT), nSamples))
names(geneTraitSignificance_E2WT) = paste("GS.", "E2WT", sep="")  # Asigna nombre adecuado
names(GSPvalue_E2WT) = paste("p.GS.", "E2WT", sep="")

geneTraitSignificance_E3WT = as.data.frame(cor(datExpr, datTraits$E3WT, use = "p"))
GSPvalue_E3WT = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_E3WT), nSamples))
names(geneTraitSignificance_E3WT) = paste("GS.", "E3WT", sep="")  # Asigna nombre adecuado
names(GSPvalue_E3WT) = paste("p.GS.", "E3WT", sep="")

geneTraitSignificance_E4WT = as.data.frame(cor(datExpr, datTraits$E4WT, use = "p"))
GSPvalue_E4WT = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_E4WT), nSamples))
names(geneTraitSignificance_E4WT) = paste("GS.", "E4WT", sep="")  # Asigna nombre adecuado
names(GSPvalue_E4WT) = paste("p.GS.", "E4WT", sep="")

# Ahora añade estas columnas al data frame geneInfo0
geneInfo0 = cbind(geneInfo0, 
                  geneTraitSignificance_E2WT, GSPvalue_E2WT,
                  geneTraitSignificance_E3WT, GSPvalue_E3WT,
                  geneTraitSignificance_E4WT, GSPvalue_E4WT)

# Ordenar los módulos por su significancia para KO (también podrías hacerlo para E2WT, E3WT, E4WT)
modOrder = order(-abs(cor(MEs, KO, use = "p")));

# Añadir la información de la membresía en el módulo en el orden elegido
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Ordenar los genes en el data frame geneInfo primero por color del módulo, luego por su significancia con KO
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.KO))  # Puedes cambiar KO por E2WT, E3WT o E4WT si prefieres
geneInfo = geneInfo0[geneOrder, ]
geneInfo<-read.table("geneInfo.microglia.txt", sep="\t",header=T)
cruchaga_genes<-c("APOD","ARFIP2","ARL1","ARL2","ATE1","ATG7","BCDIN3D","BIRC2","BTG1","C14orf93","C5orf38","CABP7","CCDC24","CCL25","CDK5RAP3","CDKN2B","CEP20","CHCHD7","CLEC4G","CLUL1","COASY","COMTD1","COPS8","CRP","CRYZL1","CTF1","DCUN1D5","DHRS9","DMKN","EFCAB1","EMILIN3","FAM20A","FBLIM1","FCHSD1","FOXO1","GCLM","GGT2","HAX1","HBQ1","HES1","HSD17B8","HYKK","IRF6","ITPK1","KCNRG","KCTD2","KRT19","LCN1","LMO4","LRRN1","MAGEA3","MAP1LC3B2","MAPK6","MAPK8","MOB4","MTMR7","NFIA","OTULIN","PICK1","PLEKHO2","PSMB10","RDH16","RND1","RNF122","S100A13","SAR1A","SELENOS","SLC27A2","SNRPF","SPATC1L","SPC25","ST8SIA1","STAR","STK10","TBCA","TIMM13","TMCC3","TP53I11","TRIM54","UNG","VPS29","WNT10B","ZNF483")


# Enrichment analysis

library(org.Hs.eg.db)
selected<-select(org.Hs.eg.db, as.character(rownames(geneInfo)), "ENTREZID", "SYMBOL")
selected2<-merge(selected,geneInfo,by.x="SYMBOL",by.y="row.names")
head(selected2)
selected2[which(selected2$SYMBOL=="APOE"),]
moduleColors<-selected2$moduleColor
GOenr = GOenrichmentAnalysis(moduleColors, selected2$ENTREZID, organism = "human", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.txt", sep = "\t", quote = TRUE, row.names = FALSE)
head(tab)
library(dplyr)
tab.table <-tab %>%
    group_by(module) %>%
    slice_max(n = 5, order_by = enrichmentP)
# tab.table[,c(1,5,6,13)]

tab.table[which(grepl("black|tan|red|purple",tab.table$module)),c(1,5,6,13)]



tan<-names(datExpr)[moduleColors=="tan"]
red<-names(datExpr)[moduleColors=="red"]
black<-names(datExpr)[moduleColors=="black"]
purple<-names(datExpr)[moduleColors=="purple"]

print("TAN")
length(tan)
sort(tan)
print("RED")
length(red)
sort(red)
print("BLACK")
length(black)
sort(black)
print("PURPLE")
length(purple)
sort(purple)

black.Enrich <- gost(
  black, 
  organism = "hsapiens", 
  ordered_query = FALSE, 
  multi_query = FALSE, 
  significant = TRUE,
  exclude_iea = FALSE, 
  measure_underrepresentation = FALSE, 
  evcodes = TRUE, 
  user_threshold = 0.05, 
  correction_method = "fdr",
  domain_scope = "annotated", 
  custom_bg = NULL, 
  numeric_ns = "", 
  sources = NULL, 
  as_short_link = FALSE
)

dfblack <- as.data.frame(black.Enrich$result)


dfblack <- data.frame(apply(dfblack, 2, as.character))

if("p_value" %in% colnames(dfblack) && "source" %in% colnames(dfblack)) {
  dfblack.table <- dfblack %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)  # Ordena y selecciona los 5 mejores resultados por p_value
  
  print(dfblack.table[, c(3, 10, 11)])
} else {
  print("Las columnas 'p_value' y/o 'source' no están presentes en el dataframe.")
}

print("**black**")
pblack <- gostplot(black.Enrich, capped = FALSE, interactive = TRUE) 
pblack

write.table(dfblack, "black.Enrichments.txt", sep="\t", row.names = FALSE)


library(gprofiler2)

red.Enrich<-gost(red, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfred<-as.data.frame(red.Enrich$result)
dfred <- data.frame(apply(dfred,2,as.character))
head(dfred[,c(3,10,11)])

dfred.table <-dfred %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfred.table[,c(3,10,11)]

print("**red**")
pred <- gostplot(red.Enrich, capped = FALSE, interactive = TRUE) 
pred

write.table(dfred,"red.Enrichments.txt", sep="\t", row.names=F)

library(gprofiler2)

purple.Enrich<-gost(purple, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfpurple<-as.data.frame(purple.Enrich$result)
dfpurple <- data.frame(apply(dfpurple,2,as.character))
head(dfpurple[,c(3,10,11)])

dfpurple.table <-dfpurple %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfpurple.table[,c(3,10,11)]

print("**purple**")
ppurple <- gostplot(purple.Enrich, capped = FALSE, interactive = TRUE) 
ppurple

write.table(dfpurple,"purple.Enrichments.txt", sep="\t", row.names=F)

geneInfo$SYMBOL<-rownames(geneInfo)
library("STRINGdb")
string_db <- STRINGdb$new( version="11.5", species=10090,
 score_threshold=200, network_type="full", input_directory="")

tan_mapped <- string_db$map(geneInfo[which(rownames(geneInfo) %in% tan),], "SYMBOL", removeUnmappedRows = TRUE )
hits.tan <- tan_mapped$STRING_id
# jpeg("tan.module.jpeg", units='cm',height=20,width=20,res=300)
# string_db$plot_network( hits.tan )
# dev.off()

tan_mapped.halo <- string_db$add_diff_exp_color( subset(tan_mapped, p.GS.KO<1),
 logFcColStr="GS.KO" )
 # post payload information to the STRING server
 payload_tan <- string_db$post_payload(tan_mapped.halo$STRING_id,
 colors=tan_mapped.halo$color )
 # display a STRING network png with the "halo"

# jpeg("tan.module.halo.jpeg", units='cm',height=20,width=20,res=300)
 string_db$plot_network( hits.tan, payload_id=payload_tan) 
# dev.off()
                     
library("STRINGdb", lib = "/home/julian/R/library")
string_db <- STRINGdb$new( version="11.5", species=10090,
 score_threshold=200, network_type="full", input_directory="")

black_mapped <- string_db$map(geneInfo[which(rownames(geneInfo) %in% black),], "SYMBOL", removeUnmappedRows = TRUE )
hits.black <- black_mapped$STRING_id
# jpeg("black.module.jpeg", units='cm',height=20,width=20,res=300)
# string_db$plot_network( hits.black )
# dev.off()

black_mapped.halo <- string_db$add_diff_exp_color( subset(black_mapped, p.GS.KO<1),
 logFcColStr="GS.KO" )
 # post payload information to the STRING server
 payload_black <- string_db$post_payload(black_mapped.halo$STRING_id,
 colors=black_mapped.halo$color )
 # display a STRING network png with the "halo"

# jpeg("black.module.halo.jpeg", units='cm',height=20,width=20,res=300)
 string_db$plot_network( hits.black, payload_id=payload_black) 
# dev.off()
# dev.off()

string_db <- STRINGdb$new( version="11.5", species=10090,
 score_threshold=200, network_type="full", input_directory="")

purple_mapped <- string_db$map(geneInfo[which(rownames(geneInfo) %in% purple),], "SYMBOL", removeUnmappedRows = TRUE )
hits.purple <- purple_mapped$STRING_id
# jpeg("purple.module.jpeg", units='cm',height=20,width=20,res=300)
# string_db$plot_network( hits.purple )
# dev.off()

purple_mapped.halo <- string_db$add_diff_exp_color( subset(purple_mapped, p.GS.KO<1),
 logFcColStr="GS.KO" )
 # post payload information to the STRING server
 payload_purple <- string_db$post_payload(purple_mapped.halo$STRING_id,
 colors=purple_mapped.halo$color )
 # display a STRING network png with the "halo"

# jpeg("purple.module.halo.jpeg", units='cm',height=20,width=20,res=300)
 string_db$plot_network( hits.purple, payload_id=payload_purple) 
# dev.off()

library("STRINGdb")
string_db <- STRINGdb$new( version="11.5", species=10090,
 score_threshold=200, network_type="full", input_directory="")

red_mapped <- string_db$map(geneInfo[which(rownames(geneInfo) %in% red),], "SYMBOL", removeUnmappedRows = TRUE )
hits.red <- red_mapped$STRING_id
# jpeg("red.module.jpeg", units='cm',height=20,width=20,res=300)
# string_db$plot_network( hits.red )
# dev.off()

red_mapped.halo <- string_db$add_diff_exp_color( subset(red_mapped, p.GS.KO<1),
 logFcColStr="GS.KO" )
 # post payload information to the STRING server
 payload_red <- string_db$post_payload(red_mapped.halo$STRING_id,
 colors=red_mapped.halo$color )
 # display a STRING network png with the "halo"

# jpeg("red.module.halo.jpeg", units='cm',height=20,width=20,res=300)
 string_db$plot_network( hits.red, payload_id=payload_red) 
# dev.off()



# Hubs
ADJ1=abs(cor(datExpr,use="p"))^15
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1[order(-Alldegrees1$kTotal) ,], n=20)

dim(datExpr)  
length(moduleColors) 

b<-chooseTopHubInEachModule(
   datExpr, 
   moduleColors, 
   omitColors = NA, 
   power = 10, 
   type = "unsigned")
b

b<-chooseTopHubInEachModule(
   datExpr, 
   moduleColors, 
   omitColors = NA, 
   power = 10, 
   type = "signed")
b


                     
