#extendedSequence <- c("ATCCATGGAGCAGTACCAGGTGGACTCGGAGCTTA",
#    "GCTTAAAGAGGACTGTTTGGACAACATGGAAGCGG",
#    "CAACATGGAAGCGGCTGTTTTTGGAGTGGTGGAGC",
#    "CATGGAAGCGGCTGTTTTTGGAGTGGTGGAGCCTC")

#featureWeightMatrix <- read.table(system.file("extdata", "Morenos-Mateo.csv",
#            package = "CRISPRseek"), sep =",", header = TRUE)

#calculategRNAEfficiencyCRISPRscan(extendedSequence, featureWeightMatrix)

#' @importFrom BiocGenerics cbind rbind do.call lapply


calculategRNAEfficiencyCRISPRscan <- function(extendedSequence, 
    featureWeightMatrix)
{
   	featureWeightMatrix[,1] = toupper(featureWeightMatrix[,1])
	efficiency <- featureWeightMatrix[featureWeightMatrix[,1] == "INTERCEPT", 2]
	features<- featureWeightMatrix[featureWeightMatrix[,1] != "INTERCEPT", 1]
	fWeights <- featureWeightMatrix[featureWeightMatrix[,1] != "INTERCEPT", 2] 
	featureNames <- gsub("[0-9]+", "", features)
	featureStart <- as.numeric(gsub("[ACGT]+", "", features))
	featureEnd <- featureStart + nchar(featureNames) - 1
	featureWeights <- cbind(featureStart, featureEnd, featureNames, fWeights)
	featureWeights <- featureWeights[order(featureWeights[,1]),]

        this.features <- numeric(dim(featureWeights)[1])
        
        all.features <- do.call(rbind, lapply(1:length(extendedSequence), 
            function(j) {
                for (i in 1:length(this.features))
                {
                    this.features[i] = as.numeric(substr(extendedSequence[j], 
                        as.numeric(featureWeights[i,1]), 
                        as.numeric(featureWeights[i,2])) == featureWeights[i,3])
                }
                this.features
        }))
	efficiency <- 100* (efficiency + all.features %*% as.numeric(featureWeights[,4]))
	efficiency
 }
