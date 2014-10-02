extendedSeq = c("TGGATTGTATAATCAGCATGGATTTGGAAC",
"TCAACGAGGATATTCTCAGGCTTCAGGTCC",
"GTTACCTGAATTTGACCTGCTCGGAGGTAA",
"CTTGGTGTGGCTTCCTTTAAGACATGGAGC",
"CATACAGGCATTGAAGAAGAATTTAGGCCT",
"AGTACTATACATTTGGCTTAGATTTGGCGG",
"TTTTCCAGATAGCCGATCTTGGTGTGGCTT",
"AAGAAGGGAACTATTCGCTGGTGATGGAGT"
)
effi.Matrix <- read.csv("/Users/zhuj/Dev/Bioconductor/Truck/CRISPRseek/inst/extdata/DoenchNBT2014.csv",header=TRUE)

calculategRNAEfficiency <- function(extendedSeq, featureWeightMatrix, gRNA.size = 20)
{
    ### add getgRNAefficiency = TRUE/FALSE to offTargetAnalysis
### If yest, need to ### add gRNAextended sequence 30bp (4bp before, and 3bp after) to offTargetAnalysis output
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
	n.C <- unlist(lapply(1:length(extendedSeq), function(i) {
						 table(factor(s2c(substr(extendedSeq[i],1,gRNA.size)), levels=c("C")))
						 }))
	n.G <- unlist(lapply(1:length(extendedSeq), function(i) {
						 table(factor(s2c(substr(extendedSeq[i],1,gRNA.size)), levels=c("G")))
						 }))
	GC.content <- n.C + n.G
	GClow.coef <- as.numeric(GC.content <10)
	GChigh.coef <- as.numeric(GC.content >10)
	efficiency <- featureWeightMatrix[featureWeightMatrix[,1] == "INTERCEPT", 2]
	efficiency <- efficiency +
		GClow.coef * (10 - GC.content) * GClow + 
		GChigh.coef * ( GC.content - 10) * GChigh
	this.features <- numeric(dim(featureWeights)[1])
	
	all.features <- do.call(rbind, lapply(1:length(extendedSeq), function(j) {
		   for (i in 1:length(this.features))
		   {
				this.features[i] = as.numeric(substr(extendedSeq[j], 
					as.numeric(featureWeights[i,1]), 
					as.numeric(featureWeights[i,2])) == featureWeights[i,3])
		   }
		   this.features
	 }))
	
	efficiency <- efficiency + all.features %*% as.numeric(featureWeights[,4])

	efficiency <- 1/(1 + exp(-efficiency))
	efficiency
}
calculategRNAEfficiency(extendedSeq, featureWeightMatrix = effi.Matrix)
