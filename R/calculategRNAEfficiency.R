calculategRNAEfficiency <- function(extendedSequence, baseBeforegRNA, featureWeightMatrix, gRNA.size = 20)
{
   	featureWeightMatrix[,1] = toupper(featureWeightMatrix[,1])
	efficiency <- featureWeightMatrix[featureWeightMatrix[,1] == "INTERCEPT", 2]
	GClow <- featureWeightMatrix[featureWeightMatrix[,1] == "GC_LOW", 2]
	GChigh <- featureWeightMatrix[featureWeightMatrix[,1] == "GC_HIGH", 2]
	features<- featureWeightMatrix[!featureWeightMatrix[,1] %in% c("INTERCEPT","GC_LOW","GC_HIGH"), 1]
	fWeights <- featureWeightMatrix[!featureWeightMatrix[,1] %in% c("INTERCEPT","GC_LOW","GC_HIGH"), 2]
	featureNames <- gsub("[0-9]+", "", features)
	featureStart <- as.numeric(gsub("[ACGT]+", "", features))
	featureEnd <- featureStart + nchar(featureNames) - 1
	featureWeights <- cbind(featureStart, featureEnd, featureNames, fWeights)
	featureWeights <- featureWeights[order(featureWeights[,1]),]
        n.cores <- detectCores() - 1
    	n.cores <- min(n.cores, length(extendedSequence))
        cl <- makeCluster(n.cores)
        clusterExport(cl, varlist = c("extendedSequence", 
            "substr", "s2c", "baseBeforegRNA", "gRNA.size", "featureWeights"),
            envir = environment()) 
	n.C <- unlist(parLapply(cl, 1:length(extendedSequence), function(i) {
	    table(factor(s2c(substr(extendedSequence[i],baseBeforegRNA + 1,
			 baseBeforegRNA + gRNA.size)), levels=c("C")))
		 }))
	n.G <- unlist(parLapply(cl, 1:length(extendedSequence), function(i) {
	    table(factor(s2c(substr(extendedSequence[i],baseBeforegRNA + 1, 
			baseBeforegRNA + gRNA.size)), levels=c("G")))
		 }))
	GC.content <- n.C + n.G
	GClow.coef <- as.numeric(GC.content <10)
	GChigh.coef <- as.numeric(GC.content >10)
	efficiency <- featureWeightMatrix[featureWeightMatrix[,1] == "INTERCEPT", 2]
	efficiency <- efficiency +
		GClow.coef * (10 - GC.content) * GClow + 
		GChigh.coef * ( GC.content - 10) * GChigh
	this.features <- numeric(dim(featureWeights)[1])
	
	all.features <- do.call(rbind, parLapply(cl, 
            1:length(extendedSequence), function(j) {
	    for (i in 1:length(this.features))
	    {
	        this.features[i] = as.numeric(substr(extendedSequence[j], 
		   as.numeric(featureWeights[i,1]), 
		   as.numeric(featureWeights[i,2])) == featureWeights[i,3])
	    }
	    this.features
	 }))
        stopCluster(cl)	
	efficiency <- efficiency + all.features %*% as.numeric(featureWeights[,4])

	efficiency <- 1/(1 + exp(-efficiency))
	efficiency
}
