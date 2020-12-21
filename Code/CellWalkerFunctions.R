#CellWalker Functions
require("TxDb.Mmusculus.UCSC.mm10.knownGene")
require("org.Mm.eg.db")

#this case the fraction of reads in marker genes scaled by expression
compute.fraction_marker = function(cellTypes, RNA_markers, ATAC_Mat){
	mmProm = promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene))
	mmGene = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
	cellType_Peak_map = sapply(cellTypes, function(x) {
		markerGenes = suppressMessages(biomaRt::select(org.Mm.eg.db, RNA_markers$Gene[RNA_markers$Cluster==x], "ENTREZID","ALIAS"))
		markerGenes = markerGenes[!is.na(markerGenes$ENTREZID) & markerGenes$ENTREZID %in% mmProm$gene_id,]
		markerPromoterLoc = mmProm[markerGenes$ENTREZID]
		markerGeneLoc = mmGene[markerGenes$ENTREZID]
		markerOverlaps = findOverlaps(as(ATAC_Peaks, "GRanges"),c(markerPromoterLoc,markerGeneLoc))
		list(peak=markerOverlaps@from,gene=markerGenes$ALIAS[(markerOverlaps@to-1)%%length(markerGenes$ALIAS)+1])
	})

	cell_peak_counts = Matrix::colSums(ATAC_Mat)
	cells_in_markers = c()
	for(cellType in cellTypes){
		peaks_in_markers = cellType_Peak_map[["peak",cellType]]
		genes_in_markers = cellType_Peak_map[["gene",cellType]]
		gene_exp = RNA_markers$`Fold Change (log)`[match(genes_in_markers,RNA_markers$Gene)] 

		cells_in_markers = cbind(cells_in_markers, Matrix::colSums(ATAC_Mat[peaks_in_markers,]*gene_exp)/cell_peak_counts)
	}

	cells_in_markers
}

#Compute jaccard on sprase matrix
jaccard = function(m) {
  A = tcrossprod(m)
  im = Matrix::which(A>0, arr.ind=TRUE)
  b = Matrix::rowSums(m!=0)
  
  Aim = A[im]
  
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  J 
}

#combine weighted label edges with cell edges
combine.graph = function(cells_in_markers, distMat, weight){
	i = dim(cells_in_markers)[2]
	rbind(cbind(matrix(0,i,i),t(weight*cells_in_markers)), cbind(weight*cells_in_markers,distMat))
}

#Solve random walk w/ restart
random.walk = function(adj, r=0.5){
  len = dim(adj)[1]

  D = Diagonal(x=Matrix::rowSums(adj))
  W = solve(D)%*%adj
  
  mat = r*solve(Diagonal(len)-r*W)
  
  mat
}

#Compute cell homogeneity
compute.cell.homogeneity = function(cellTypes, labelsAdjType, infMat){
	i = length(cellTypes)
  log(median(unlist(sapply(cellTypes, function(cellType) 
    	median(sapply(which(labelsAdjType==cellType), function(x) {
    		median(infMat[x+i,which(labelsAdjType==cellType)+i])/median(infMat[x+i,which(!labelsAdjType==cellType)+i])
    	}))
  ))))
}

#Compute uncertain labels
compute.uncertain.labels = function(cellTypes, cellTypeInf){
	#difference between top two label scores
	uncertaintyScore = apply(cellTypeInf, 1, function(x) sort(x, decreasing = TRUE)[1]-sort(x, decreasing = TRUE)[2])
	uncertaintyTypes = apply(cellTypeInf, 1, function(x) c(cellTypes[order(x, decreasing = TRUE)][1], cellTypes[order(x, decreasing = TRUE)][2]))
	#all pairs of labels where difference in scores is in the bottom 10%
	uncertainPairs = uncertaintyTypes[,uncertaintyScore < quantile(uncertaintyScore, .1)]
	uncertainMat = sapply(cellTypes, function(c1) sapply(cellTypes, function(c2) length(which((uncertainPairs[1,]==c1 & uncertainPairs[2,]==c2) | (uncertainPairs[1,]==c2 & uncertainPairs[2,]==c1)))))
	uncertainMat
}

#Compare two cell types
compare.types = function(type1, type2, cellLabels, cellTypes, cellTypeInf){
  cellWalk$normMat[cellLabels==type1 | cellLabels==type2,which(cellTypes==type1)]-cellWalk$normMat[cellLabels==type1 | cellLabels==type2,which(cellTypes==type2)]
}

#Correlate cell types with TADS
correlate.TADS = function(TADRanges, type1, type2, scores, cellLabels, ATAC_Mat, ATAC_Peaks){
  PeakRanges = as(ATAC_Peaks,"GRanges")
  normSum = colSums(ATAC_Mat[,cellLabels==type1 | cellLabels==type2])
  TADCor = sapply(1:length(TADRanges), function(TAD) {
    whichPeaks = which(countOverlaps(PeakRanges, TADRanges[TAD])>0)
    if(length(whichPeaks)==0){NA}
    else if(length(whichPeaks)==1){cor.test(ATAC_Mat[whichPeaks,cellLabels==type1 | cellLabels==type2]/normSum, scores)$estimate}
    else{cor.test(colSums(ATAC_Mat[whichPeaks,cellLabels==type1 | cellLabels==type2])/normSum, scores)$estimate}
  })
  TADCor
}


#Map bulk peaks
labelBulk = function(bulkPeaks, infMat, ATAC_Mat, ATAC_Peaks, cellTypes){
  
  ATAC_Peaks = as(ATAC_Peaks, "GRanges")
  peakOverlaps = findOverlaps(ATAC_Peaks, bulkPeaks)
  
  infCellOnType = sapply(1:length(bulkPeaks), function(e) {
      whichPeaks = peakOverlaps@from[peakOverlaps@to==e]
      if(length(whichPeaks)==0){
        rep(NA, length(cellTypes))
      } else{
        testCells = ATAC_Mat[whichPeaks,]
        if(length(whichPeaks)==1){
          whichCells = which(testCells>0)
        } else{
          whichCells = which(colSums(testCells)>0)
        }
        rowSums(infMat[1:length(cellTypes),length(cellTypes)+whichCells])/length(whichCells)
      }
  })
  
  mappedLabel = apply(infCellOnType, 2, function(x) ifelse(length(which(is.na(x)))==0,cellTypes[order(x, decreasing = TRUE)[1]],NA))
  mappedLabel
}
