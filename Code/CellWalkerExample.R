#CellWalker Example:
#Run CellWalker on scATAC-seq portion of SNARE-seq data

require("data.table")
require("Matrix")
source("CellWalkerFunctions.R")

#Load scATAC-seq data, GEO accession: GSE126074
ATAC_Barcodes = fread("GSE126074_AdBrainCortex_SNAREseq_chromatin.barcodes.tsv.gz", header=FALSE)$V1
ATAC_Peaks = fread("GSE126074_AdBrainCortex_SNAREseq_chromatin.peaks.tsv.gz", header=FALSE)$V1
ATAC_Mat = readMM("GSE126074_AdBrainCortex_SNAREseq_chromatin.counts.mtx.gz")
rownames(ATAC_Mat) = ATAC_Peaks
colnames(ATAC_Mat) = ATAC_Barcodes
ATAC_Mat = ATAC_Mat>0 #binarize

#Load marker genes from "High-throughput sequencing of the transcriptome and chromatin accessibility in the same cell," Supplementary Table 1
RNA_markers = fread("adult_cerebral_cortex_markers.csv")
cellTypes = names(table(RNA_markers$Cluster))

#Generate a label to cell mapping (label-cell edge weights)
cells_in_markers = compute.fraction_marker(cellTypes, RNA_markers, ATAC_Mat)

#compute jaccard distance between all cells (cell-cell edge weights)
distMat = jaccard(t(ATAC_Mat))

#diffuse labels over graph for different label edge weights
i = length(cellTypes)
cellHomogeneityList = list()
cellTypeInfList = list()
typeTypeInfList = list()
cellLabelList = list()
for(weight in 10^seq(-2,4,1)){
	#combined graph
	simMat = combine.graph(cells_in_markers, distMat, weight)
		
	#compute influence with random walk
	infMat = random.walk(simMat)

	#Label-to-cell influence
	cellTypeInf = as.matrix(infMat[-(1:i),(1:i)])
	colnames(cellTypeInf) = cellTypes

	#Label-to-label influence
	typeTypeInf = as.matrix(infMat[(1:i),(1:i)])
	colnames(typeTypeInf) = cellTypes
	rownames(typeTypeInf) = cellTypes

	#normalize label-to-cell scores for ranking
    labelsAdj = apply(cellTypeInf, 2, function(x) (x-mean(x))/sd(x))
    labelsAdj = apply(labelsAdj, 2, function(x) x/max(x))
    labelsAdjType = apply(labelsAdj, 1, function(x) cellTypes[order(x, decreasing = TRUE)][1])
    
    #compute cell homogeneity and store portions of information matrix
    infMat = as.matrix(infMat)
    cellHomogeneityList[[as.character(weight)]] =  compute.cell.homogeneity(cellTypes, labelsAdjType, infMat)
    cellTypeInfList[[as.character(weight)]] = labelsAdj
    typeTypeInfList[[as.character(weight)]] = typeTypeInf
    cellLabelList[[as.character(weight)]] = labelsAdjType
}

#plot cell homogenity across edge weights
require("ggplot2")
cellHomogeneityScores = unlist(cellHomogeneityList)
p = ggplot() + geom_line(aes(10^seq(-2,4,1), cellHomogeneityScores)) + scale_x_log10() + scale_y_log10() + theme_classic() + ylab("Predicted Cell Homogeneity") + xlab("Label Edge Weight")
ggplot(p)

#select optimal edge weight
weight = as.character(10^seq(-2,4,1)[order(cellHomogeneityScores, decreasing=TRUE)[1]])

#heirachical clustering of cell types based on label-label weights
typeTypeInf = typeTypeInfList[[weight]]
plot(hclust(dist(typeTypeInf)))

#look at uncertain labels
cellTypeInf = cellTypeInfList[[weight]]
uncertainMat = compute.uncertain.labels(cellTypes, cellTypeInf)
uncertainMatMelt = melt(uncertainMat)
uncertainMatMelt$Var1 = factor(uncertainMatMelt$Var1, levels = cellTypes)
uncertainMatMelt$Var2 = factor(uncertainMatMelt$Var2, levels = cellTypes)
ggplot(uncertainMatMelt, aes(Var1, Var2)) + geom_tile(aes(fill = value), colour = "black") +scale_fill_gradient(low = "white",high = "steelblue") + theme(axis.text.x = element_text(angle = 45, hjust=0), axis.title.x=element_blank(), axis.text.y = element_text(angle = 45), axis.title.y=element_blank()) + scale_x_discrete(position = "top") 

#most likely cell labels
cellLabels = cellLabelList[[weight]]

#compare two very similar cell types
l34Score = compare.types("Ex-L3/4-Rorb", "Ex-L3/4-Rmst", cellLabels, cellTypes, cellTypeInf)
ggplot() + geom_density(aes(l34Score)) + xlab("Ex-L3/4 Rorb vs Rmst Score")

#compare two less similar cell types
l34_23Score = compare.types("Ex-L3/4-Rorb", "Ex-L2/3-Rasgrf2", cellLabels, cellTypes, cellTypeInf)
ggplot() + geom_density(aes(l34_23Score)) + xlab("Ex-L3/4 Rorb vs Ex-L2/3 Rasgrf2 Score")

#analyze TADs
TADRanges = fread("total.combined.domain")
TADRanges = GRanges(TADRanges$V1, IRanges(TADRanges$V2,TADRanges$V3))
ATAC_Mat = as.matrix(ATAC_Mat) #unsparsify matrix for faster calculations
l34_TADCor = correlate.TADS(TADRanges, "Ex-L3/4-Rorb", "Ex-L3/4-Rmst", l34Score, cellLabels, ATAC_Mat, ATAC_Peaks)
ggplot() + geom_density(aes(l34_TADCor)) + xlab("Ex-L3/4 Rorb vs Rmst Score - Accessibility Correlation in TADs")

#map enhancer
Vista_forebrain = fread("ENCFF225JUL.bed.gz")
Vista_forebrain = GRanges(Vista_forebrain$V1, IRanges(Vista_forebrain$V2,Vista_forebrain$V3))
enhancer_labels = labelBulk(Vista_forebrain, infMat, ATAC_Mat, ATAC_Peaks, cellType)
enhancer_labels = enhancer_labels[!is.na(enhancer_labels)] # some can't me mapped due to no overlap between bulk and single cell data
ggplot() + geom_bar(aes(enhancer_labels)) + xlab("Cell Type") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))




