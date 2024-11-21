# Helper functions used by offTargetAnalysis()
# Validate input parameters passed to offTargetAnalysis()
# Stop and exit if any error detected in input parameters
#' @importFrom rlang inform warn
#' @importFrom openxlsx read.xlsx write.xlsx
#' @importFrom Biostrings readDNAStringSet
#' @author Kai Hu

# Read in potential.gRNAs from path or DNAStringSet:
readgRNAsFromFile <- function(inputFilePath = NULL, format = NULL) {
  if (class(inputFilePath) != "DNAStringSet") {
    if (!file.exists(inputFilePath)) {
      stop("inputfile specified as ", inputFilePath, " does not exists!")
    }
    if (format == "fasta" || format == "fastq") {
      potential.gRNAs <- readDNAStringSet(inputFilePath, format, use.names = TRUE)
    } else {
      stop("format needs to be either fasta,fastq or bed!")
    }
  } else {
    potential.gRNAs <- inputFilePath
    if (length(names(potential.gRNAs)) == 0) {
      names(potential.gRNAs) <- paste("gRNAs", 1:length(potential.gRNAs), sep = "")
    }
  }
  return(potential.gRNAs)
}

# Return a set of grouped argument lists:
getArgGroups <- function(inputFilePath = NULL,
                         overlap.gRNA.positions = NULL,
                         baseEditing = NULL, 
                         targetBase = NULL, 
                         editingWindow = NULL,
                         editingWindow.offtargets = NULL,
                         primeEditing = NULL,
                         findPairedgRNAOnly = NULL,
                         annotatePaired = NULL,
                         paired.orientation = NULL,
                         pairOutputFile = NULL,
                         PAM = NULL,
                         PAM.location = NULL,
                         PAM.size = NULL,
                         gRNA.pattern = NULL,
                         gRNA.size = NULL,
                         min.gap = NULL,
                         max.gap = NULL,
                         name.prefix = NULL,
                         format = NULL,
                         rule.set = NULL,
                         chrom_acc = NULL,
                         PBS.length = NULL,
                         RT.template.length = NULL,
                         RT.template.pattern = NULL,
                         targeted.seq.length.change = NULL,
                         bp.after.target.end = NULL,
                         target.start = NULL,
                         target.end = NULL,
                         primeEditingPaired.output = NULL,
                         corrected.seq = NULL,
                         calculategRNAEfficacy = TRUE,
                         featureWeightMatrixFile = NULL,
                         baseBeforegRNA = NULL,
                         baseAfterPAM = NULL,  
                         efficacyFile = NULL,
                         enable.multicore = NULL,
                         n.cores.max = NULL,
                         BSgenomeName = NULL, #
                         genomeSeqFile = NULL,
                         txdb = NULL,
                         chromToSearch = NULL, #
                         chromToExclude = NULL, #
                         max.mismatch = NULL, #
                         PAM.pattern = NULL, #
                         allowed.mismatch.PAM = NULL, #
                         outputDir = NULL, ##
                         exportAllgRNAs = NULL, ##
                         gRNAoutputName = NULL, ##
                         findgRNAs = NULL, ###### add paired and RE info to summary, get potential gRNA and filter gRNA can not use it
                         findgRNAsWithREcutOnly = NULL, ### filtergRNAs
                         REpatternFile = NULL, ###
                         minREpatternSize = NULL, ###
                         subPAM.position = NULL, #### getOfftargetScoreWrap 
                         subPAM.activity = NULL, ####
                         scoring.method = NULL, ####
                         mismatch.activity.file = NULL, ####
                         weights = NULL, ####
                         fetchSequence = NULL, ##### getOfftargetSummary
                         orgAnn = NULL, #####
                         ignore.strand = NULL, #####
                         min.score = NULL, #####
                         topN = NULL, #####
                         topN.OfftargetTotalScore = NULL, #####
                         upstream = NULL, #####
                         downstream = NULL, #####
                         annotateExon = NULL, #####
                         calculategRNAefficacyForOfftargets = NULL, #####
                         DNA_bulge = NULL, ###### getOfftargetWithBulge
                         RNA_bulge = NULL, ######
                         method.findOffTargetsWithBulge = NULL) {
  findgRNAs_args_core <- list(inputFilePath = inputFilePath,
                              overlap.gRNA.positions = overlap.gRNA.positions,
                              baseEditing = baseEditing, 
                              targetBase = targetBase, 
                              editingWindow = editingWindow,
                              primeEditing = primeEditing,
                              findPairedgRNAOnly = findPairedgRNAOnly,
                              annotatePaired = annotatePaired,
                              paired.orientation = paired.orientation,
                              pairOutputFile = pairOutputFile,
                              PAM = PAM,
                              PAM.location = PAM.location,
                              PAM.size = PAM.size,
                              gRNA.pattern = gRNA.pattern,
                              gRNA.size = gRNA.size,
                              min.gap = min.gap,
                              max.gap = max.gap,
                              name.prefix = name.prefix, #*
                              format = format,
                              rule.set = rule.set,
                              chrom_acc = chrom_acc)
  findgRNAs_args_prime <- list(PBS.length = PBS.length,
                               RT.template.length = RT.template.length,
                               RT.template.pattern = RT.template.pattern,
                               targeted.seq.length.change = targeted.seq.length.change,
                               bp.after.target.end = bp.after.target.end,
                               target.start = target.start,
                               target.end = target.end,
                               primeEditingPaired.output = primeEditingPaired.output,
                               corrected.seq = corrected.seq)
  findgRNAs_args_efficacy <- list(calculategRNAEfficacy = TRUE,
                                  featureWeightMatrixFile = featureWeightMatrixFile,
                                  baseBeforegRNA = baseBeforegRNA,
                                  baseAfterPAM = baseAfterPAM,  
                                  efficacyFile = efficacyFile,
                                  enable.multicore = enable.multicore,
                                  n.cores.max = n.cores.max)
  searchHits_args <- list(chromToSearch = chromToSearch, #
                          chromToExclude = chromToExclude, #
                          max.mismatch = max.mismatch, #
                          PAM.pattern = PAM.pattern, #
                          allowed.mismatch.PAM = allowed.mismatch.PAM, #
                          PAM = PAM, 
                          PAM.size = PAM.size,
                          PAM.location = PAM.location,
                          gRNA.size = gRNA.size, 
                          baseEditing = baseEditing, 
                          targetBase = targetBase,
                          editingWindow = editingWindow.offtargets)
  outputPotentialgRNAs_args <- list(outputDir = outputDir, 
                                    exportAllgRNAs = exportAllgRNAs, 
                                    gRNAoutputName = gRNAoutputName, 
                                    inputFilePath = inputFilePath, 
                                    format = format, 
                                    overlap.gRNA.positions = overlap.gRNA.positions)
  filtergRNAs_args <- list(findgRNAsWithREcutOnly = findgRNAsWithREcutOnly,
                           REpatternFile = REpatternFile,
                           format = format,
                           minREpatternSize = minREpatternSize,
                           overlap.gRNA.positions = overlap.gRNA.positions)
  getgRNASummary_args <- list(PAM = PAM, 
                              PAM.size = PAM.size,
                              PAM.location = PAM.location, 
                              gRNA.size = gRNA.size, 
                              outputDir = outputDir)
  getOfftargetScoreWrap_args <- list(gRNA.size = gRNA.size,
                                     canonical.PAM = PAM,
                                     PAM.size = PAM.size,
                                     PAM.location = PAM.location,
                                     subPAM.position = subPAM.position,
                                     subPAM.activity = subPAM.activity,
                                     scoring.method = scoring.method,
                                     mismatch.activity.file = mismatch.activity.file,
                                     weights = weights)
  getOfftargetSummary_args <- list(outputDir = outputDir, 
                                   BSgenomeName = BSgenomeName, 
                                   genomeSeqFile = genomeSeqFile,
                                   txdb = txdb,
                                   gRNA.size = gRNA.size,
                                   PAM.location = PAM.location, 
                                   PAM.size = PAM.size,
                                   featureWeightMatrixFile = featureWeightMatrixFile,
                                   rule.set = rule.set, 
                                   chrom_acc = chrom_acc,
                                   fetchSequence = fetchSequence,
                                   orgAnn = orgAnn, 
                                   ignore.strand = ignore.strand,
                                   min.score = min.score, 
                                   topN = topN,
                                   topN.OfftargetTotalScore = topN.OfftargetTotalScore,
                                   upstream = upstream, 
                                   downstream = downstream,
                                   annotateExon = annotateExon,
                                   baseBeforegRNA = baseBeforegRNA,
                                   baseAfterPAM = baseAfterPAM, 
                                   calculategRNAefficacyForOfftargets = calculategRNAefficacyForOfftargets)
  getOfftargetWithBulge_args <- list(PAM.size = PAM.size,
                                     PAM.pattern = PAM.pattern,
                                     PAM.location = PAM.location,
                                     max.mismatch = max.mismatch,
                                     DNA_bulge = DNA_bulge, 
                                     RNA_bulge = RNA_bulge,
                                     BSgenomeName = BSgenomeName, 
                                     genomeSeqFile = genomeSeqFile,
                                     chromToSearch = chromToSearch,
                                     cas_offinder_version = method.findOffTargetsWithBulge)
  filterCasOffinder_args = list(max.mismatch = max.mismatch, 
                                gRNA.size = gRNA.size,
                                PAM = PAM,
                                PAM.size = PAM.size,
                                PAM.pattern = PAM.pattern, 
                                PAM.location = PAM.location,
                                allowed.mismatch.PAM = allowed.mismatch.PAM, 
                                baseEditing = baseEditing,
                                targetBase = targetBase,
                                editingWindow = editingWindow)
  
  return(list(findgRNAs_args_core, 
              findgRNAs_args_prime, 
              findgRNAs_args_efficacy, 
              searchHits_args, 
              outputPotentialgRNAs_args,
              filtergRNAs_args,
              getgRNASummary_args,
              getOfftargetScoreWrap_args,
              getOfftargetSummary_args,
              getOfftargetWithBulge_args,
              filterCasOffinder_args))
}

# Get all potential gRNAs:
getPotentialgRNAs <- function(findgRNAs = NULL, 
                              chromToSearch = NULL, 
                              inputFilePath = NULL,
                              format = NULL,
                              useEfficacyFromInputSeq = NULL, 
                              findgRNAs_args_core = NULL, 
                              findgRNAs_args_prime = NULL, 
                              findgRNAs_args_efficacy = NULL,
                              outputPotentialgRNAs_args = NULL) {
  inform("Searching for gRNAs ...")
  
  if (findgRNAs) {
    if (is.null(chromToSearch) || useEfficacyFromInputSeq) {
      potential.gRNAs <- do.call(CRISPRseek::findgRNAs, c(findgRNAs_args_core, findgRNAs_args_prime, findgRNAs_args_efficacy))
    } else {
      potential.gRNAs <- do.call(CRISPRseek::findgRNAs, c(findgRNAs_args_core, findgRNAs_args_prime))
    }
  } else {
    potential.gRNAs <- readgRNAsFromFile(inputFilePath, format = format)
  }
  
  if (findgRNAs) {
    do.call(outputPotentialgRNAs, c(list(potential.gRNAs), outputPotentialgRNAs_args))
  }
  inform("Searching for gRNAs: done!")
  return(potential.gRNAs)
}

# Output potentila gRNAs to file:
outputPotentialgRNAs <- function(potential.gRNAs = NULL, 
                                 outputDir = NULL, 
                                 exportAllgRNAs = NULL, 
                                 gRNAoutputName = NULL, 
                                 inputFilePath = NULL, 
                                 format = NULL, 
                                 overlap.gRNA.positions = NULL) {
  if (length(potential.gRNAs) > 0) {
    if ((exportAllgRNAs == "fasta" || exportAllgRNAs == "all")) {
      writeXStringSet(potential.gRNAs, 
                      filepath = file.path(outputDir,
                                           paste(gRNAoutputName, "allgRNAs.fa", sep = "")))
    }
    
    if ((exportAllgRNAs == "genbank" || exportAllgRNAs == "all")) {
      if (class(inputFilePath) == "DNAStringSet") {
        subjects <- inputFilePath
      } else {
        subjects <- readDNAStringSet(inputFilePath, format=format, use.names = TRUE)
      }

      names(subjects) <- gsub( "\t", "", names(subjects))
      names(subjects) <- gsub( "\n", "", names(subjects))
      names(subjects) <- gsub( " ", "", names(subjects))
      locuses <- names(subjects)
      
      names.gRNA <- names(potential.gRNAs)
      for (i in 1:length(locuses)) {
        thisLocus <- gsub("'", "", locuses[i])
        thisLocus <- gsub(" ", "", thisLocus)
        thisSeq <- tolower(as.character(subjects[[i]]))
        n.bp <- nchar(thisSeq)
        temp <- strsplit(names.gRNA, split = paste(thisLocus, "_gR", sep = ""))
        locus <- paste("LOCUS       ", thisLocus,
                       "                     ", n.bp,
                       " bp    dna     linear   UNK", sep = "")
        definition <- paste("DEFINITION  CRISPRseek output for ",
                            gRNAoutputName, " sequence", sep = "")
        accession <- "ACCESSION   unknown"
        features <- "FEATURES             Location/Qualifiers"
        header = rbind(locus, definition, accession, features)
        found.gRNA <- 0
        for (j in 1:length(temp)) {
          if (length(temp[[j]]) > 1) {
            found.gRNA <- found.gRNA + 1
            if (found.gRNA == 1) {
              thisFile <- file.path(outputDir, paste(thisLocus, "gbk", sep = "."))
              write(header, thisFile)
            }
            if (length(grep("f", temp[[j]])) > 0) {
              temp1 <-strsplit(temp[[j]], "f")
              isForward <- TRUE
            } else {
              temp1 <-strsplit(temp[[j]], "r")
              isForward <- FALSE
            }
            feature <- temp1[[2]][2]
            feature[is.na(feature)] <- ""
            location <- temp1[[2]][1]
            if (isForward) {
              Start <- location
              End <- as.numeric(Start) + max(overlap.gRNA.positions) - min(overlap.gRNA.positions)
              write(paste("     misc_bind       ", Start, "..",
                          End, sep = ""), append = TRUE, sep = "\n",
                    file = thisFile)
              write(paste("                     /note=\"gRNAf",
                          as.character(feature),
                          "\"", sep = ""), 
                    append = TRUE, sep="\n", 
                    file = thisFile)
            } else {
              End <- location
              Start <- as.numeric(End) - max(overlap.gRNA.positions) + min(overlap.gRNA.positions)
              write(paste("     misc_bind       complement(",
                          Start, "..", End, ")", sep = ""), 
                    append = TRUE, sep="\n", file = thisFile)
              write(paste("                     /note=\"gRNAr",
                          feature,	"\"", sep = ""), 
                    append = TRUE, sep = "\n", file = thisFile)
            }
          }
        }
        if (found.gRNA > 0) {
          write("ORIGIN", append = TRUE, sep = "\n", file = thisFile)
          seq.lines <- floor(nchar(thisSeq) / 60) + 1
          for (k in 1:seq.lines) {
            line.start <- (k - 1) * 60 + 1
            line.end <- min(line.start + 59, nchar(thisSeq))
            n.leading.spaces <- 9 - nchar(line.start)
            leading.spaces <- paste(rep(" ", n.leading.spaces), collapse = "")
            seq.thisLine <- substr(thisSeq, line.start, line.end)
            len.thisLine <- nchar(seq.thisLine)
            n.seg <- floor(len.thisLine / 10) + 1
            for (l in 1:n.seg) {
              seg.start <- (l -1) * 10 + 1
              seg.end <- min(seg.start + 9, len.thisLine)
              if (l == 1) {
                seq.thisLine.formatted <- substr(seq.thisLine, seg.start, seg.end)
              } else {
                seq.thisLine.formatted <- paste( seq.thisLine.formatted,
                                                 substr(seq.thisLine, seg.start, seg.end),
                                                 sep = " ")
              }
            }
            write(paste(leading.spaces, line.start, " ",
                        seq.thisLine.formatted, sep = ""),
                  append = TRUE, sep="\n", file = thisFile)
          }
          write("//", append = TRUE, sep="\n", file = thisFile)
        }
      }
    }
  }
}

# Get and output filtered gRNA to file:
getFiltergRNAs <- function(potential.gRNAs = NULL,
                           findgRNAs = NULL,
                           filtergRNAs = NULL, 
                           filtergRNAs_args = NULL, 
                           findgRNAsWithREcutOnly = NULL, 
                           annotatePaired = NULL,
                           findPairedgRNAOnly = FALSE, 
                           pairOutputFile = NULL,
                           REcutDetailFile = NULL) {
  if (length(potential.gRNAs) == 0) {
    return(c("", ""))
  } else {
    if (findPairedgRNAOnly) {
      gRNAs.RE <- do.call(filtergRNAs, c(list(potential.gRNAs), list(pairOutputFile = pairOutputFile), filtergRNAs_args))
      REcutDetails  <- gRNAs.RE$gRNAREcutDetails
      REcutDetails$gap <- as.numeric(REcutDetails$gap)
      if ("REcutStart" %in% names(REcutDetails)) {
        REcutDetails$REcutStart <- as.numeric(REcutDetails$REcutStart)
        REcutDetails$REcutEnd <- as.numeric(REcutDetails$REcutEnd)
      }
      REcutDetails$ForwardREcutStart <- as.numeric(REcutDetails$ForwardREcutStart)
      REcutDetails$ForwardREcutEnd <- as.numeric(REcutDetails$ForwardREcutEnd)
      REcutDetails$ReverseREcutStart <- as.numeric(REcutDetails$ReverseREcutStart)
      REcutDetails$ReverseREcutEnd <- as.numeric(REcutDetails$ReverseREcutEnd)
      write.xlsx(REcutDetails[order(as.character(REcutDetails$ForwardgRNAName)), ], file = REcutDetailFile, rowNames = FALSE)
    } else {
      gRNAs.RE <- do.call(filtergRNAs, c(list(potential.gRNAs), list(pairOutputFile = ""), filtergRNAs_args))
      REcutDetails  <- gRNAs.RE$gRNAREcutDetails
      if ("gap" %in% names(REcutDetails)) {
        REcutDetails$gap <- as.numeric(REcutDetails$gap)
      }
      if ("REcutStart" %in% names(REcutDetails)) {
        REcutDetails$REcutStart <- as.numeric(REcutDetails$REcutStart)
        REcutDetails$REcutEnd <- as.numeric(REcutDetails$REcutEnd)
      }
      write.xlsx(REcutDetails[order(as.character(REcutDetails$REcutgRNAName)), ], file = REcutDetailFile, rowNames = FALSE)
    }

    if (findgRNAsWithREcutOnly) {
      gRNAs  <- gRNAs.RE$gRNAs.withRE
    } else {
      gRNAs <- potential.gRNAs
    }
    
    pairedInformation <- ""
    if (findgRNAs) {
      if (annotatePaired || findPairedgRNAOnly) {
        pairedInformation <- read.xlsx(pairOutputFile)
      }
    } 
    
    names(gRNAs) <- gsub("[\t\n ]", "", names(gRNAs))
    return(list(gRNAs, pairedInformation, REcutDetails))
  }
}

# Get and output gRNA summary (if no off targets detected)
getgRNASummary <- function(gRNAs = NULL, 
                           PAM = NULL, 
                           PAM.size = NULL, 
                           PAM.location = NULL, 
                           gRNA.size = NULL, 
                           outputDir = NULL) {
  if (PAM.location == "3prime") {
    x <- paste(substr(as.character(gRNAs), 1, gRNA.size), PAM, sep ="")
  } else {
    x <- paste(PAM, substr(as.character(gRNAs), PAM.size + 1, gRNA.size + PAM.size), sep ="")
  }
  
  summary <- cbind(names = names(gRNAs), 
                   gRNAsPlusPAM = x,
                   top5OfftargetTotalScore = rep("NA", length(gRNAs)),
                   top10OfftargetTotalScore = rep("NA", length(gRNAs)),
                   top1Hit.onTarget.MMdistance2PAM = rep("perfect match not found", length(gRNAs)))
  write.xlsx(summary, file = file.path(outputDir, "Summary.xlsx"), rowNames = FALSE)
  return(summary)
}

# Get offtarget scores
getOfftargetScoreWrap <- function(hits = NULL, 
                                  gRNA.size = NULL,
                                  canonical.PAM = NULL,
                                  PAM.size = NULL,
                                  PAM.location = NULL,
                                  subPAM.position = NULL, 
                                  subPAM.activity = NULL, 
                                  scoring.method = NULL, 
                                  mismatch.activity.file = NULL,
                                  weights = NULL) {
  inform("Building feature vectors for scoring ...")
  featureVectors <- buildFeatureVectorForScoring(hits = hits,
                                                 gRNA.size = gRNA.size,
                                                 canonical.PAM = canonical.PAM,
                                                 subPAM.position = subPAM.position,
                                                 PAM.location = PAM.location, 
                                                 PAM.size = PAM.size)
  inform("Building feature vectors built for scoring: done!")
  
  inform("Calculating off-target scores ...")
  if (scoring.method ==  "CFDscore") {
    scores <- getOfftargetScore2(featureVectors,
                                 subPAM.activity = subPAM.activity,
                                 mismatch.activity.file = mismatch.activity.file)
  } else {
    scores <- getOfftargetScore(featureVectors, weights = weights)
  }
  inform("Calculating off-target scores: done!")
  return(scores)
}

# Get offtarget scores: for off-target with bulges
getOfftargetScoreBulgeWrap <- function(hits = NULL, 
                                       gRNA.size = NULL,
                                       canonical.PAM = NULL,
                                       PAM.size = NULL,
                                       PAM.location = NULL,
                                       subPAM.position = NULL, 
                                       subPAM.activity = NULL, 
                                       scoring.method = NULL, 
                                       mismatch.activity.file = NULL,
                                       weights = NULL) {
  if (dim(hits)[1] > 0) {
    inform("Building feature vectors (bulge) for scoring ...")
    # alns <- convertHits2Alns(hits = hits)
    hits <- addAlnInfoToHits(hits = hits, PAM.size = PAM.size, PAM.location = PAM.location, gRNA.size = gRNA.size)
    featureVectors <- buildFeatureVectorForScoringBulge(alns = hits,
                                                        gRNA.size = gRNA.size,
                                                        canonical.PAM = canonical.PAM,
                                                        PAM.location = PAM.location, 
                                                        PAM.size = PAM.size,
                                                        insertion.symbol = "^")
    inform("Building feature vectors (bulge) for scoring: done!")
    inform("Calculating off-target scores ...")
    if (scoring.method ==  "CFDscore") {
      scores <- getOfftargetScore2(featureVectors,
                                   subPAM.activity = subPAM.activity,
                                   mismatch.activity.file = mismatch.activity.file)
    } else {
      scores <- getOfftargetScore(featureVectors, weights = weights)
    }
    inform("Calculating off-target scores: done!")
  } else {
    inform("Skipping off-target with bulge!")
    dimnames = list(c(), c("gRNAPlusPAM_bulge", "OffTargetSequence_bulge", "n.RNABulge", "n.DNABulge", "gRNA.insertion", "gRNA.deletion", "offTarget_sequence", "PAM.sequence", "pos.mismatch", "pos.insertion", "pos.deletion", paste0("IsInsertion.pos", 1:20), paste0("IsDeletion.pos", 1:20), paste0("IsMismatch.pos", 1:20), "strand", "chrom", "chromStart", "chromEnd", "name", "n.mismatch", "forViewInUCSC", "alignment", "NGG", "mean.neighbor.distance.mismatch", "mismatch.type", "subPAM", "score"))
    scores <- data.frame(matrix(ncol = length(dimnames[[2]]), nrow = 0, dimnames = dimnames))
  }

  return(scores)
}

# Annotate and filter offtargets 
getOfftargetSummary <- function(scores = NULL,
                                outputDir = NULL, 
                                BSgenomeName = NULL, 
                                genomeSeqFile = NULL,
                                fetchSequence = NULL,
                                txdb = NULL,
                                orgAnn = NULL, 
                                ignore.strand = NULL,
                                min.score = NULL, 
                                topN = NULL,
                                topN.OfftargetTotalScore = NULL,
                                upstream = NULL, 
                                downstream = NULL,
                                annotateExon = NULL,
                                baseBeforegRNA = NULL,
                                baseAfterPAM = NULL, 
                                gRNA.size = NULL,
                                PAM.location = NULL, 
                                PAM.size = NULL,
                                featureWeightMatrixFile = NULL,
                                rule.set = NULL, 
                                chrom_acc = NULL,
                                calculategRNAefficacyForOfftargets = NULL) {
  offTargets <- filterOffTarget(scores = scores, 
                                outputDir = outputDir,
                                BSgenomeName = BSgenomeName, 
                                genomeSeqFile = genomeSeqFile,
                                fetchSequence = fetchSequence,
                                txdb = txdb,
                                orgAnn = orgAnn, 
                                ignore.strand = ignore.strand,
                                min.score = min.score, 
                                topN = topN,
                                topN.OfftargetTotalScore = topN.OfftargetTotalScore,
                                upstream = upstream, 
                                downstream = downstream,
                                annotateExon = annotateExon, 
                                baseBeforegRNA = baseBeforegRNA,
                                baseAfterPAM = baseAfterPAM, 
                                gRNA.size = gRNA.size,
                                PAM.location = PAM.location,
                                PAM.size = PAM.size,
                                featureWeightMatrixFile = featureWeightMatrixFile,
                                rule.set = rule.set, 
                                chrom_acc = chrom_acc,
                                calculategRNAefficacyForOfftargets = calculategRNAefficacyForOfftargets)
  inform("Annotating and filtering off-targets ...")
  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"))
  if (dim(summary)[2] == 1) {
    summary <- as.data.frame(t(data.matrix(offTargets$summary)))
  }
    
  inform("Annotating and filtering off-targets: done!")
  return(list(offTargets, summary))
}

# Add paired info and RE to the summary
addInfoSummary <- function(summary = NULL,
                           findgRNAs = NULL, 
                           annotatePaired = NULL,
                           findPairedgRNAOnly = NULL,
                           pairedInformation = NULL,
                           REcutDetails = NULL) {
  if (findgRNAs && (annotatePaired || findPairedgRNAOnly)) {
    inform("Adding paired information ...")
    PairedgRNAName <- unlist(lapply(1:dim(summary)[1], function(i) {
      as.character(gsub("^\\s+|\\s+$", "",
                        paste(unique(pairedInformation[as.character(pairedInformation$ForwardgRNAName) == as.character(summary$names[i]),]$ReversegRNAName),
                              unique(pairedInformation[as.character(pairedInformation$ReversegRNAName) == as.character(summary$names[i]),]$ForwardgRNAName), collapse = " ")))
    }))
    inform("Adding paired information: done!")
  }

  if (findPairedgRNAOnly && findgRNAs) {
    inform("Adding RE information ...")
    REname <- unlist(lapply(1:dim(summary)[1], function(i) {
      gsub("^\\s+|\\s+$", "", 
           gsub("NA", "", paste(unique(REcutDetails[as.character(REcutDetails$ForwardREcutgRNAName) == as.character(summary$names[i]),]$ForwardREname),
                                unique(REcutDetails[as.character(REcutDetails$ReverseREcutgRNAName) == as.character(summary$names[i]), ]$ReverseREname), collapse = " ")))
      }))
    summary <- cbind(summary, PairedgRNAName, REname)
    inform("Adding RE information: done!")
  } else {
    REname <- unlist(lapply(1:dim(summary)[1], function(i) {
      gsub("^\\s+|\\s+$", "", gsub("NA", "", paste(unique(REcutDetails[as.character(REcutDetails$REcutgRNAName) == as.character(summary$names[i]), ]$REname), collapse = " ")))
    }))
    summary <- cbind(summary, REname)
  }
  return(summary)
}

# Write gRNAs to bed
writegRNAsToBed <- function(summary = NULL,
                            offTargets = NULL,
                            bedFile = NULL,
                            useScore = NULL,
                            overlap.gRNA.positions = NULL) {
  inform("Writting gRNAs to bed file ...")
  
  on.target <- offTargets$offtargets
  on.target <- unique(subset(on.target, on.target$n.mismatch == 0 & on.target$isCanonicalPAM == 1))
  gRNA.bed <- ""
  if (dim(on.target)[1] == 0) {
    warn("No on-target found for the input gRNAs with your search criteria!")
  } else {
    gRNA.bed <- unique(cbind(as.character(on.target$chrom),
                             as.character(on.target$chromStart),
                             as.character(on.target$chromEnd), 
                             as.character(on.target$name),	
                             as.numeric(as.character(on.target$gRNAefficacy)) * 1000,
                             as.character(on.target$strand),
                             as.character(on.target$chromStart),
                             as.character(on.target$chromEnd)))
    if (!useScore) {
      gRNA.bed <- cbind(gRNA.bed, rep("255,0,0", dim(gRNA.bed)[1]))
      gRNA.bed[gRNA.bed[,6] == "-",9] = "0,255,0"
    }

    gRNA.bed[, 2] = as.numeric(gRNA.bed[, 2]) - 1 # UCSC genome browser uses 0-based coordinates
    gRNA.bed[, 3] = as.numeric(gRNA.bed[, 3])
    gRNA.bed[gRNA.bed[,6] == "+" ,7] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "+" ,2]) + min(overlap.gRNA.positions) - 1
    gRNA.bed[gRNA.bed[,6] == "-" ,7] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "-" ,3]) - max(overlap.gRNA.positions)
    gRNA.bed[gRNA.bed[,6] == "+", 8] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "+", 2]) + max(overlap.gRNA.positions)
    gRNA.bed[gRNA.bed[,6] == "-", 8] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "-", 3]) - min(overlap.gRNA.positions) + 1
    
    write.table("track name=\"gRNA sites\" 
	               description=\"CRISPRseek\" visibility=2 useScore=1 itemRgb=\"On\"", 
                file = bedFile, col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(gRNA.bed, file = bedFile, sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
    inform("Writting gRNAs to bed file: done!")
  }
  return(gRNA.bed)
}


# Scan for REsites in flanking region
scanREsites <- function(summary = NULL,
                        offTargets = NULL,
                        BSgenomeName = NULL,
                        REcutDetails = NULL,
                        outputUniqueREs = NULL) {
  inform("Scanning REsites in flanking regions ...")
  on.target <- offTargets$offtargets
  on.target <- unique(subset(on.target, on.target$n.mismatch == 0 & on.target$isCanonicalPAM == 1))
  
  if (dim(on.target)[1] > 0 && outputUniqueREs && !is.null(BSgenomeName) && class(BSgenomeName) == "BSgenome") {
    summary$REs.isUnique100 <- uniqueREs(REcutDetails = REcutDetails, summary = summary, offTargets$offtargets, scanUpstream = 100, scanDownstream = 100, BSgenomeName = BSgenomeName)
    summary$REs.isUnique50 <- uniqueREs(REcutDetails = REcutDetails, summary = summary, offTargets$offtargets, scanUpstream = 50, scanDownstream = 50, BSgenomeName = BSgenomeName)
  } else {
    summary$REs.isUnique100 <- ""
    summary$REs.isUnique50 <- ""
  }
  summary$uniqREin200 <- summary$REs.isUnique100 
  summary$uniqREin100 <- summary$REs.isUnique50 
  inform("Scanning REsites in flanking regions: done!")
  return(summary)
}

# Update summary: add on.targets, RE
updateSummary <- function(summary = NULL, 
                          outputDir = NULL, 
                          offTargets = NULL,
                          useEfficacyFromInputSeq = NULL,
                          efficacyFile = NULL) {
  on.target <- getOntarget(offTargets)
  # on.target <- offTargets$offtargets
  # on.target <- unique(subset(on.target, on.target$n.mismatch == 0 & on.target$isCanonicalPAM == 1))
  
  if (dim(on.target)[1] == 0) {
    warn("No on-target found for the input gRNAs with your search criteria!")
  } else {
    on.target <- unique(cbind(as.character(on.target$name),
                              as.character(on.target$forViewInUCSC),
                              as.character(on.target$extendedSequence),
                              as.numeric(on.target$gRNAefficacy)))
    colnames(on.target) = c("names", "forViewInUCSC", "extendedSequence", "gRNAefficacy")
    
    if (useEfficacyFromInputSeq) {
      on.target <- as.data.frame(on.target[, 1:2])
      inputEfficacy <- read.table(efficacyFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      inputEfficacy <- as.data.frame(cbind(name = inputEfficacy$name, extendedSequence = inputEfficacy$extendedSequence, gRNAefficacy = inputEfficacy$gRNAefficacy))
      on.target <- merge(on.target, inputEfficacy, by.x = "names", by.y ="name")
    }
    
    summary <- unique(merge(on.target, summary, by = "names", all = TRUE))
    summary$gRNAefficacy <- as.numeric(summary$gRNAefficacy)
    summary$uniqREin200 <- rep("", dim(summary)[1])
    summary$uniqREin100 <- rep("", dim(summary)[1])
    # summary <- cbind(summary, uniqREin200 = REs.isUnique100, uniqREin100 = REs.isUnique50)
    # summary$uniqREin200 <- as.character(summary$uniqREin200)
    # summary$uniqREin100 <- as.character(summary$uniqREin100)
  }
  write.xlsx(summary[order(as.character(summary$names)), ], file = file.path(outputDir, "Summary.xlsx"), rowNames = FALSE)
  return(summary)
}

# Update summary: add foldgRNA
updateSummaryFoldgRNAs <- function(foldgRNAs = NULL, 
                                   summary = NULL, 
                                   gRNA.size = NULL, 
                                   gRNA.backbone = NULL, 
                                   temperature = NULL) {
  if (foldgRNAs) {
    source(system.file("extdata/foldgRNAs.R", package = "CRISPRseek"))
    gRNAs.withoutPAM <- substr(as.character(summary$gRNAsPlusPAM), 1, gRNA.size)
    folded.gRNAs <- foldgRNAs(gRNAs.withoutPAM, gRNA.backbone = gRNA.backbone,
                              temperature = temperature)
    if (length(dim(folded.gRNAs)) > 0) {
      if (dim(folded.gRNAs)[1] >1) {
        summary <- cbind(summary, folded.gRNAs[,-1])
      } else {
        summary <- data.frame(c(summary, folded.gRNAs[,-1]))
      }
    }
  }
  return(summary)
}

# Update summary: add no.target gRNAs
updateSummaryNoTargetgRNAs <- function(summary = NULL,
                                       gRNAs = NULL, 
                                       gRNA.size = NULL,
                                       PAM = NULL,
                                       PAM.size = NULL,
                                       PAM.location = NULL,
                                       on.target = NULL,
                                       outputDir = NULL) {
  ### even there is no perfect target for a gRNA, it will be kept in the summary file
  ### need to calculate the topN offtarget score and distance correctly yet if include those gRNAs without target
  
  gRNAs.notInGenome <- setdiff(names(gRNAs), summary$names)
  if (length(gRNAs.notInGenome) > 0) {
    dat2 <- data.frame(matrix(nrow = length(gRNAs.notInGenome), ncol = dim(summary)[2]))
    colnames(dat2) <- colnames(summary)
    dat2$names <- gRNAs.notInGenome
    
    if (PAM.location == "3prime") {
      dat2$gRNAsPlusPAM <- paste0(substr(as.character(gRNAs[names(gRNAs) %in% gRNAs.notInGenome]), 1, gRNA.size), PAM)
    } else {
      dat2$gRNAsPlusPAM <- paste0(PAM, substr(as.character(gRNAs[names(gRNAs) %in%gRNAs.notInGenome]), PAM.size + 1, + PAM.size + gRNA.size))
    }

    dat2$top1Hit.onTarget.MMdistance2PAM <- rep("perfect match not found", length(gRNAs.notInGenome))
    summary <- rbind(summary, dat2)
  }
  
  if (dim(on.target)[1] == 0) {
    write.xlsx(summary[order(as.character(summary$names)), ], file = file.path(outputDir, "Summary.xlsx"), rowNames = FALSE)
  } else {
    write.xlsx(summary[order(as.character(summary$forViewInUCSC)), ], file = file.path(outputDir, "Summary.xlsx"), rowNames = FALSE)
  }

  return(summary)
}

# Update summary: add Lindel results
addLindelRes <- function(predIndelFreq = NULL,
                         predictIndelFreq.onTargetOnly = NULL,
                         offTargets = NULL,
                         on.target = NULL,
                         BSgenomeName = NULL,
                         baseBeforegRNA.indelFreq = NULL,
                         baseAfterPAM.indelFreq = NULL,
                         method.indelFreq = NULL,
                         summary = NULL,
                         outputDir = NULL) {
  indelFreq <- ""
  fs2 <- ""
  
  if (predIndelFreq) {
    if (predictIndelFreq.onTargetOnly) {
      targets <- unique(subset(offTargets$offtargets, offTargets$offtargets$n.mismatch == 0 & offTargets$offtargets$isCanonicalPAM == 1))
    } else {
      targets <- subset(offTargets$offtargets, offTargets$offtargets$isCanonicalPAM == 1)
    }
    
    extendedSequence <- getExtendedSequence(targets, BSgenomeName = BSgenomeName, baseBeforegRNA =  baseBeforegRNA.indelFreq, baseAfterPAM = baseAfterPAM.indelFreq, forMethod = method.indelFreq)

    tryCatch((indelFreqFS <- predictRelativeFreqIndels(extendedSequence, 
                                                       method = method.indelFreq)), 
             error = function(e) {print(e); })
    
    indelFreqFS <- tryCatch(
      predictRelativeFreqIndels(extendedSequence, method = method.indelFreq),
      error = function(e) {
        message(e$message)
        return(NULL)
      }
    )
    
    
    if (!is.null(indelFreqFS)) {
      fs <- unlist(lapply(indelFreqFS, function(x) { x$fs }))
      indelFreq <- lapply(indelFreqFS, function(x) { x$indel })
      names(indelFreq) <- paste(targets[,1], targets[,2], targets[,3], sep = ",")
      
      entropy <- unlist(lapply(indelFreq, function(x) {
        if (length(x) > 1) {
          sum(-as.numeric(x[,2])/100 * log(as.numeric(x[,2])/100, base = 450), na.rm = TRUE)
        } else {
          NA
        }
      }))
      
      fs2 <- data.frame(cbind(names = as.character(targets[,1]), frameshift = fs, entropy = entropy, n.mismatch = as.character(targets$n.mismatch)))
      fs2[, 1] <- as.character(fs2[, 1])
      summary <- data.frame(summary)
      summary[, 1] <- as.character(summary[, 1])
      summary <- merge(subset(fs2, fs2[, 4] == 0)[, -4], summary, all.y = TRUE)
      
      write.xlsx(summary[order(as.character(summary$forViewInUCSC)), ], file = file.path(outputDir, "Summary.xlsx"), rowNames = FALSE)
      
      if (!predictIndelFreq.onTargetOnly) {
        offTargets$offtargets[,3] <- as.character(offTargets$offtargets[,3])
        fs3 <- cbind(OffTargetSequence =  as.character(targets[,3]), frameshift = fs, entropy = entropy)
        targets <- merge(offTargets$offtargets, fs3, all.x = TRUE)
        offTargets$offtargets <- targets
        write.xlsx(targets, file = file.path(outputDir, "OfftargetAnalysis.xlsx"), rowNames = FALSE)
      }
    }
  }
  
  return(list(summary, indelFreq, fs2))
}

# Combine scores and scores_bulge:
combineScoresBulge <- function(scores = NULL, scores_bulge = NULL) {
  # combine scores_bulge with scores 
  names(scores_bulge)[names(scores_bulge) == "offTarget_sequence"] <- "OffTargetSequence"
  scores_bulge$gRNAPlusPAM <- scores$gRNAPlusPAM[match(scores_bulge$name, scores$name)]
  scores_bulge <- scores_bulge[, !names(scores_bulge) %in% "PAM.sequence"]
  scores_bulge$mismatch.distance2PAM <- rep("", dim(scores_bulge)[1])
  
  scores$gRNAPlusPAM_bulge <- rep("", dim(scores)[1])
  scores$OffTargetSequence_bulge <- rep("", dim(scores)[1])
  scores$gRNA.insertion <- rep("", dim(scores)[1])
  scores$gRNA.deletion <- rep("", dim(scores)[1])
  scores$pos.insertion <- rep("", dim(scores)[1])
  scores$pos.deletion <- rep("", dim(scores)[1])
  scores$pos.mismatch <- apply(scores, 1, function(row) {
    res <- c()
    for (n in names(row)) {
      if (startsWith(n, "IsMismatch.pos")) {
        if (as.integer(row[n])) {
          res <- c(res, as.integer(row[n]))
        }
      }
    }
    if (is.null(res)) {
      res <- character()
    }
    return(list(res))
  })
  for (col in setdiff(names(scores_bulge), names(scores))) {
    scores[[col]] <- rep(0, dim(scores)[1])
  }

  scores_bulge <- scores_bulge[, colnames(scores)]
  scores <- rbind(scores, scores_bulge)
  
  # Remove duplicate on.target hits that are from bulged search:
  perfect_match <- grepl("^\\.+$", scores$alignment)
  with_bulge <- nchar(scores$gRNAPlusPAM_bulge) > 0
  scores <- scores[!(perfect_match & with_bulge), ]

  return(scores)
}

# Search Hits using custom reference fasta file:
searchHitsFa <- function(gRNAs = NULL,
                         gRNA.size = NULL,
                         genomeSeqFile = NULL,
                         chromToSearch = NULL,
                         chromToExclude = NULL, 
                         max.mismatch = NULL, 
                         PAM = NULL,
                         PAM.pattern = NULL,
                         PAM.size = NULL,
                         PAM.location = NULL,
                         allowed.mismatch.PAM = allowed.mismatch.PAM,
                         baseEditing = NULL, 
                         targetBase = NULL,
                         editingWindow = NULL) {
  genomeSeq <- readDNAStringSet(genomeSeqFile)
  
  if (length(chromToSearch) == 1 && tolower(chromToSearch) == "all") {
    chromInd <- 1:length(genomeSeq)
  } else {
    chromInd <-  which(names(genomeSeq) %in% chromToSearch)
  }
  chromInd <-  setdiff(chromInd, which(names(genomeSeq) %in% chromToExclude))
  
  for (j in chromInd) {
    if (j == chromInd[1]) {
      hits <- searchHits(gRNAs = gRNAs,
                         PAM = PAM, 
                         PAM.pattern = PAM.pattern,
                         seqs = genomeSeq[[j]],
                         seqname = names(genomeSeq)[j],
                         max.mismatch = max.mismatch, 
                         PAM.size = PAM.size,
                         gRNA.size = gRNA.size, 
                         allowed.mismatch.PAM = allowed.mismatch.PAM,
                         PAM.location = PAM.location,
                         baseEditing = baseEditing, 
                         targetBase = targetBase,
                         editingWindow = editingWindow,
                         outfile = tempfile(tmpdir = getwd()))
    } else {
      hits <- rbind(hits, searchHits(gRNAs = gRNAs,
                                     PAM = PAM, 
                                     PAM.pattern = PAM.pattern,
                                     seqs = genomeSeq[[j]], 
                                     seqname = names(genomeSeq)[j],
                                     max.mismatch = max.mismatch, 
                                     PAM.size = PAM.size,
                                     gRNA.size = gRNA.size, 
                                     allowed.mismatch.PAM = allowed.mismatch.PAM,
                                     PAM.location = PAM.location,
                                     baseEditing = baseEditing, 
                                     targetBase = targetBase,
                                     editingWindow = editingWindow,
                                     outfile = tempfile(tmpdir = getwd())))
    }
  }
  
  return(hits)
}

# Filter out and keep only the on.target from offTargets search result
getOntarget <- function(offTargets = NULL) {
  on.target <- offTargets$offtargets
  on.target <- on.target[grepl("^\\.+$", on.target$alignment), ]
  on.target <- unique(subset(on.target, on.target$n.mismatch == 0 & on.target$isCanonicalPAM == 1))
}
