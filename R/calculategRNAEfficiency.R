#' Calculate gRNA Efficiency
#' 
#' Calculate gRNA Efficiency for a given set of sequences and feature weight
#' matrix
#' 
#' 
#' @param extendedSequence Sequences containing gRNA plus PAM plus flanking
#' sequences. Each sequence should be long enough for building features
#' specified in the featureWeightMatrix
#' @param baseBeforegRNA Number of bases before gRNA used for calculating gRNA
#' efficiency, default 4
#' @param featureWeightMatrix a data frame with the first column containing
#' significant features and the second column containing the weight of
#' corresponding features. In the following example, DoenchNBT2014 weight
#' matrix is used. Briefly, features include 
#' \itemize{
#' \item{INTERCEPT}
#' \item{GC_LOW} - {penalty for low
#' GC content in the gRNA sequence}
#' \item{GC_HIGH} - {penalty for high GC content in
#' the gRNA sequence}
#' \item{G02} - {means G at second position of the
#' extendedSequence}
#' \item{GT02} - {means GT di-nucleotides starting at 2nd position of
#' the extendedSequence}
#' }
#' To understand how is the feature weight matrix is
#' identified, or how to use alternative feature weight matrix file, please see
#' Doench et al., 2014 for details.
#' @param gRNA.size The size of the gRNA, default 20
#' @param enable.multicore Indicate whether enable parallel processing, default
#' FALSE. For super long sequences with lots of gRNAs, suggest set it to TRUE
#' @param n.cores.max Indicating maximum number of cores to use in multi core
#' mode, i.e., parallel processing, default 6. Please set it to 1 to disable
#' multicore processing for small dataset.
#' @return DNAStringSet consists of potential gRNAs that can be input to
#' filtergRNAs function directly
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references Doench JG, Hartenian E, Graham DB, Tothova Z, Hegde M, Smith I,
#' Sullender M, Ebert BL, Xavier RJ, Root DE. Rational design of highly active
#' sgRNAs for CRISPR-Cas9-mediated gene inactivation. Nat Biotechnol. 2014 Sep
#' 3. doi: 10.1038 nbt.3026
#' http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
#' @keywords misc
#' @examples
#' 
#'  	extendedSequence <- c("TGGATTGTATAATCAGCATGGATTTGGAAC",
#' 		"TCAACGAGGATATTCTCAGGCTTCAGGTCC",
#' 		"GTTACCTGAATTTGACCTGCTCGGAGGTAA",
#' 		"CTTGGTGTGGCTTCCTTTAAGACATGGAGC",
#' 		"CATACAGGCATTGAAGAAGAATTTAGGCCT",
#' 		"AGTACTATACATTTGGCTTAGATTTGGCGG",
#' 		"TTTTCCAGATAGCCGATCTTGGTGTGGCTT",
#' 		"AAGAAGGGAACTATTCGCTGGTGATGGAGT"
#' 	)
#' 	featureWeightMatrixFile <- system.file("extdata", "DoenchNBT2014.csv", 
#' 		package = "CRISPRseek")
#' 	featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
#' 	calculategRNAEfficiency(extendedSequence, baseBeforegRNA = 4, 
#' 		featureWeightMatrix, gRNA.size = 20)
#' @importFrom BiocGenerics cbind order unlist table do.call rbind lapply
#' @importFrom parallel detectCores makeCluster clusterExport parLapply
#' @importFrom seqinr s2c
#' @export
calculategRNAEfficiency <- function(extendedSequence, baseBeforegRNA, 
    featureWeightMatrix, gRNA.size = 20, enable.multicore = FALSE,
    n.cores.max = 6)
{
   	featureWeightMatrix[,1] = toupper(featureWeightMatrix[,1])
	efficiency <- featureWeightMatrix[featureWeightMatrix[,1] == "INTERCEPT", 2]
	GClow <- featureWeightMatrix[featureWeightMatrix[,1] == "GC_LOW", 2]
	GChigh <- featureWeightMatrix[featureWeightMatrix[,1] == "GC_HIGH", 2]
	features<- featureWeightMatrix[!featureWeightMatrix[,1] %in% 
        c("INTERCEPT","GC_LOW","GC_HIGH"), 1]
	fWeights <- featureWeightMatrix[!featureWeightMatrix[,1] %in% 
        c("INTERCEPT","GC_LOW","GC_HIGH"), 2]
	featureNames <- gsub("[0-9]+", "", features)
	featureStart <- as.numeric(gsub("[ACGT]+", "", features))
	featureEnd <- featureStart + nchar(featureNames) - 1
	featureWeights <- cbind(featureStart, featureEnd, featureNames, fWeights)
	featureWeights <- featureWeights[order(featureWeights[,1]),]
    n.cores <- detectCores() - 1
    n.cores <- min(n.cores, length(extendedSequence))
    n.cores <- min(n.cores, n.cores.max)  
    if (enable.multicore && n.cores > 1)
    {
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
        
        all.features <- do.call(rbind, parLapply(cl, 1:length(extendedSequence), 
            function(j) {
                for (i in 1:length(this.features))
                {
                    this.features[i] = as.numeric(substr(extendedSequence[j], 
                    as.numeric(featureWeights[i,1]), 
                    as.numeric(featureWeights[i,2])) == featureWeights[i,3])
                }
                this.features
        }))
    } #### if enable.multicore is TRUE
   	else
   	{
       n.C <- unlist(lapply(1:length(extendedSequence), function(i) {
   	        table(factor(s2c(substr(extendedSequence[i],baseBeforegRNA + 1,
   	            baseBeforegRNA + gRNA.size)), levels=c("C")))
   	    }))
   	    n.G <- unlist(lapply(1:length(extendedSequence), function(i) {
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
   	} #### if enable.multicore is FALSE
	efficiency <- efficiency + all.features %*% as.numeric(featureWeights[,4])

	efficiency <- 1/(1 + exp(-efficiency))
	efficiency
}
