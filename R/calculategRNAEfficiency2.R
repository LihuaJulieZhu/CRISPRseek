#Sys.setenv(PATH = paste("/Users/ZHUJ/anaconda2/bin", Sys.getenv("PATH"), sep=":")) 
#system("python --version")
#Python 2.7.15 :: Anaconda, Inc.

#extendedSequence="TCTGGTGTCCTCCACACCAGAATCAGGGGT"

#CRISPRseek:::calculategRNAEfficiency2(extendedSequence,aa.cut = -1, per.peptide = -1)

#' @importFrom reticulate py_discover_config
#' @importFrom BiocGenerics lapply

calculategRNAEfficiency2 <- function(extendedSequence,  
     aa.cut = -1, per.peptide = -1)
{
    system2("python", args = "--version", stderr="pythonVersion.txt")
    #if(grep( "2.7", read.table("pythonVersion.txt", sep="", header=FALSE)[1,2]))
    if(grep( "2.7", py_discover_config()$version)) 
    {	
 #	unlink("pythonVersion.txt")
	origDir <- getwd()
	pythonDir <- system.file("extdata/Rule_Set_2_scoring_v1/analysis/",package = "CRISPRseek")
        origDir <- getwd()
        if (dir.exists(pythonDir))
        {
            setwd(pythonDir)
	    tryCatch((
	        efficiency <- unlist(lapply(extendedSequence, function(thisSeq) {
                   system(paste("python rs2_score_calculator.py -seq ", thisSeq, " -aa-cut ", aa.cut, " -per-peptide ", per.peptide, sep=""), intern=TRUE)
                  }))), error = function(e) {print(e); })
            efficiency <- unlist(lapply(efficiency, function(temp) {strsplit(temp, ": ")[[1]][2]}))
	    setwd(origDir)
        }
        else
        {
             stop("Rule_Set_2_scoring_v1 python code not found in the extdata directory of CRISPRseek, please update CRISPRseek to the newest version first!") 
        }
   }
   else
   {
       stop("Python 2.7 is required for calculating gRNA efficacy using Rule_Set_2_scoring published in JG. Doench, et al., Nbt, Jan 2016!")
   }
   efficiency
}
