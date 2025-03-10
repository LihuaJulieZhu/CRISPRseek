# Input validations
# Validate input parameters passed to offTargetAnalysis()
# Stop and exit if any error detected in input parameters
#' @importFrom rlang inform abort
# @author Kai Hu

fixPAMpattern <- function(PAM.pattern, PAM.location) {
  PAM.p.letters <- strsplit(PAM.pattern, split = "")[[1]]
  first_letter <- PAM.p.letters[1]
  last_letter <- PAM.p.letters[length(PAM.p.letters)]
  if (PAM.location == "3prime") {
    if (first_letter == "^") {
      abort("PAM.pattern does not match with PAM.location!")
    } else if (last_letter != "$") {
      PAM.pattern <- paste0(PAM.pattern, "$")
    }
  } else if (PAM.location == "5prime") {
    if (last_letter == "$") {
      abort("PAM.pattern does not match with PAM.location!")
    } else if (first_letter != "^") {
      PAM.pattern <- paste0("^", PAM.pattern)
    }
  }
  
  return(PAM.pattern)
}

checkPAMsize <- function(PAM, PAM.size) {
  # width_pam_pattern <- PAM.pattern %>% str_remove("^\\^") %>% str_remove("\\$$") %>% width() # PAM.pattern can contain other special characters like |, [].
  if (width(PAM) != PAM.size) {
    abort("PAM length does not match with PAM.size!")
  }
}

fixRuleSet <- function(rule.set, scoring.method, subPAM.activity, PAM.location, PAM.size, gRNA.size) {
  if (rule.set == "DeepCpf1") {
    baseBeforegRNA <- 8
    baseAfterPAM <- 26
    if (scoring.method == "CFDscore" && subPAM.activity$TT < 1) {
      subPAM.activity = hash(AA = 0,
                             AC = 0,
                             AG = 0,
                             AT = 0.1,
                             CA = 0,
                             CC = 0,
                             CG = 0,
                             CT = 0.05,
                             GA = 0,
                             GC = 0,
                             GG = 0,
                             GT = 0.05,
                             TA = 0.2,
                             TC = 0.1,
                             TG = 0.1,
                             TT = 1)
    }
  } else if (rule.set %in% c("Root_RuleSet1_2014", "Root_RuleSet2_2016", "CRISPRscan")) {
    if (PAM.location == "3prime") {
      baseBeforegRNA <- 4
      baseAfterPAM <- 3
    } else {
      baseBeforegRNA <- 4 + PAM.size
      baseAfterPAM <- 3 + gRNA.size
    }
  }
  
  return(list(baseBeforegRNA, baseAfterPAM, subPAM.activity))
}

validateScoringMethod <- function(scoring.method, mismatch.activity.file, weights, gRNA.size) {
  if (scoring.method ==  "CFDscore") {
    mismatch.activity <- read.csv(mismatch.activity.file)
    required.col <- c("Mismatch.Type", "Position", "Percent.Active")
    if (length(intersect(colnames(mismatch.activity), required.col)) != length(required.col)) {
      abort("Please rename the mismatch activity file column to contain at least these three column names: Mismatch.Type, Position, Percent.Active")
      }
  } else if (scoring.method == "Hsu-Zhang") {
    if (length(weights) !=  gRNA.size)
    abort("Please make sure the size of weights vector equals to the gRNA.size!")
    }
}

checkREpatternFile <- function(findgRNAsWithREcutOnly, findgRNAs, REpatternFile) {
  if (findgRNAsWithREcutOnly && findgRNAs && !file.exists(REpatternFile)) {
    abort("Please specify an REpattern file as fasta file with restriction enzyme recognition sequences!")
  }
}

checkDependency <- function(foldgRNAs = NULL) {
  if (foldgRNAs) {
    inform("Please ensure `GeneRfold` package is installed!")
    # comment out to suppress the warning message during R CMD check
    # tryCatch(GeneRfold::fold("TTAGCTCAATTGGTAAAGACCCTAGGCGAAGCTTAGAGGTCGCCGGTT"), 
    #          error = function(e) abort("'foldgRNAs' is set to TRUE, but package 'GeneRfold' is not available!"))
  }
}

prepOutputDir <- function(outputDir, overwrite) {
  if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != .Platform$file.sep) {
    outputDir <- paste(outputDir, "", sep = .Platform$file.sep)
  }
  if (!file.exists(outputDir)) {
    dir.create(outputDir)
  } else if (!overwrite) {
    inform(paste0(outputDir, " already exists. Please type 1 if you want to overwrite the outputDir and 2 if you want to exit."))
    
    input <- readline()
    if (input != 1) { 
      stop("Please change the outputDir!")
    }
  }
}

prepOtherFiles <- function(annotatePaired, findPairedgRNAOnly, inputFilePath, outputDir, gRNAoutputName) {
  pairOutputFile <- NULL
  if (annotatePaired || findPairedgRNAOnly) {
    pairOutputFile <- file.path(outputDir, "pairedgRNAs.xlsx")
  }
  REcutDetailFile <- file.path(outputDir, "REcutDetails.xlsx")
  bedFile <- file.path(outputDir, "gRNAsCRISPRseek.bed")
  
  if (is.null(gRNAoutputName)) {
    if (class(inputFilePath) == "DNAStringSet") {
      # below seems not necessary, indeed, the gRNAoutputName parameter seems redundant
      # assign a default "test" to gRNAoutputName
      abort("Please enter a name for the gRNA output file ('gRNAoutputName') when a DNAStringSet object is provided as 'inputFilePath'!")
    } else {
      gRNAoutputName = strsplit(basename(inputFilePath), split = ".", fixed = TRUE)[[1]][1]
    }
  }
  
  return(list(pairOutputFile, REcutDetailFile, bedFile, gRNAoutputName))
}
