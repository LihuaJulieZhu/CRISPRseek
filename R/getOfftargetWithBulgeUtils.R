#' @noRd
#' @title Obtain off-targets with bulge utility functions
#' @importFrom reticulate py_config py_discover_config py_available install_miniconda source_python
#' @importFrom gtools permutations
#' @importFrom Biostrings readDNAStringSet
#' @importFrom rlang abort
#' @importFrom stringr str_count str_sub str_to_upper str_to_lower str_replace_all str_detect str_split str_c
#' @author Kai Hu

# Helper functions for `getOfftargetWithBulge()`.

# Convert gRNA into DNAStringSet given offTargetAnalysis output list
getgRNAFromoffTargetAnalysis <- function(res) {
  gRNA_PAM <- DNAStringSet(res$summary$gRNAsPlusPAM)
  names(gRNA_PAM) <- res$summary$names
  return(gRNA_PAM)
}

# Translate regex PAM pattern to IUPAC pattern
translatePAMpattern <- function(PAM, direction = c("toIUPAC", "fromIUPAC")) {
  if (is.null(PAM)) {
    PAM <- ""
  }
  PAM <- toupper(PAM)
  if (direction == "fromIUPAC") {
    PAM <- gsub("W", "[A|T]", PAM)
    PAM <- gsub("M", "[A|C]", PAM)
    PAM <- gsub("R", "[A|G]", PAM)
    PAM <- gsub("Y", "[C|T]", PAM)
    PAM <- gsub("S", "[G|C]", PAM)
    PAM <- gsub("D", "[A|G|T]", PAM)
    PAM <- gsub("H", "[A|C|T]", PAM)
    PAM <- gsub("V", "[A|C|G]", PAM)
    PAM <- gsub("K", "[T|U|G]", PAM)
    PAM <- gsub("B", "[C|G|T]", PAM)
    PAM <- gsub("N", "[A|C|T|G]", PAM)
  } else if (direction == "toIUPAC") {
    regex_all <- lapply(2:4, function(i) {
      res <- permutations(4, i, v = c("A", "T", "C", "G"))
      apply(res, 1, function(row) {
        paste0("[", paste(row, collapse = "|"), "]")
      }
      )
    })
    for (v in regex_all) {
      if (width(v[[1]]) == 5) { # contains two nucleotides: e.g. [A|T]
        lapply(v, function(pattern) {
          if (grepl("A", pattern) && grepl("T", pattern)) {
            PAM <<- gsub(pattern, "W", PAM, fixed = TRUE)
          } else if (grepl("A", pattern) && grepl("C", pattern)) {
            PAM <<- gsub(pattern, "M", PAM, fixed = TRUE)
          } else if (grepl("A", pattern) && grepl("G", pattern)) {
            PAM <<- gsub(pattern, "R", PAM, fixed = TRUE)
          } else if (grepl("C", pattern) && grepl("T", pattern)) {
            PAM <<- gsub(pattern, "Y", PAM, fixed = TRUE)
          } else if (grepl("G", pattern) && grepl("C", pattern)) {
            PAM <<- gsub(pattern, "S", PAM, fixed = TRUE)
          }
        })
      } else if (width(v[[1]]) == 7) {
        lapply(v, function(pattern) {
          if (grepl("A", pattern) && grepl("G", pattern) && grepl("T", pattern)) {
            PAM <<- gsub(pattern, "D", PAM, fixed = TRUE)
          } else if (grepl("A", pattern) && grepl("C", pattern) && grepl("T", pattern)) {
            PAM <<- gsub(pattern, "H", PAM, fixed = TRUE)
          } else if (grepl("A", pattern) && grepl("C", pattern) && grepl("G", pattern)) {
            PAM <<- gsub(pattern, "V", PAM, fixed = TRUE)
          } else if (grepl("T", pattern) && grepl("U", pattern) && grepl("G", pattern)) {
            PAM <<- gsub(pattern, "K", PAM, fixed = TRUE)
          } else if (grepl("C", pattern) && grepl("G", pattern) && grepl("T", pattern)) {
            PAM <<- gsub(pattern, "B", PAM, fixed = TRUE)
          }
        })
      } else if (width(v[[1]]) == 9) {
        lapply(v, function(pattern) {
          PAM <<- gsub(pattern, "N", PAM, fixed = TRUE)
        })
      }
    }
  } else {
    abort("Wrong direction supplied! Choose from c('toIUPAC', 'fromIUPAC')")
  }
  return(PAM)
}

# Check and download CasOFFinder to bin/
getCasOFFinder <- function(version = c("2.4.1", "3.0.0b3")) {
  version <- match.arg(version)
  os <- Sys.info()["sysname"]
  arch <- ifelse(.Machine$sizeof.pointer == 8, "64", "32")
  if (arch != "64") {
    abort("Only 64-bit systems are supported!")
  }
  base_url <- "https://github.com/snugel/cas-offinder/releases/download/"
  base_path <- file.path(system.file(package = "CRISPRseek"), "bin", "cas-offinder")
  
  if (version == "3.0.0b3") {
    suffix <- "_64.zip"
    suffix_mac <- "_mac_x86_64.zip"
    subdir <- "build"
  } else {
    suffix <- "-64.zip"
    suffix_mac <- "_macos_x86-64.zip"
    subdir <- ""
  }
  
  if (os == "Darwin") {
    zip_url <- paste0(base_url, version, "/cas-offinder", suffix_mac)
    zip_path <- file.path(base_path, version, paste0("cas-offinder", suffix_mac))
    binary_path <- file.path(base_path, version, subdir, "cas-offinder")
  } else if (os == "Linux") {
    zip_url <- paste0(base_url, version, "/cas-offinder_linux_x86", suffix)
    zip_path <- file.path(base_path, version, paste0("cas-offinder_linux_x86", suffix))
    binary_path <- file.path(base_path, version, subdir, "cas-offinder")
  } else if (os == "Windows") {
    zip_url <- paste0(base_url, version, "/cas-offinder_windows_x86", suffix)
    zip_path <- file.path(base_path, version, paste0("cas-offinder_windows_x86", suffix))
    binary_path <- file.path(base_path, version, subdir, "cas-offinder.exe")
  } else {
    abort(paste0("Unsupported operating system: ", os))
  }
  
  if (!file.exists(binary_path)) {
    dir.create(dirname(binary_path), recursive = TRUE, showWarnings = FALSE)
    download.file(zip_url, zip_path, mode = "wb")
    unzip(zip_path, exdir = dirname(zip_path))
    system(paste0("chmod 755 ", shQuote(binary_path)))
    # inform("Cas-OFFinder binary downloaded!")
  } else {
    # inform(paste0("Cas-OFFinder already downloaded at ", binary_path))
  }
  return(binary_path)
}

# Verify CasOFFinder and determine if GPU available (returns "G" for gpu or "C" for cpu)
checkCasOFFinder <- function(binary_path) {
  status <- tryCatch({
    system2(binary_path, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    return(paste("Error:", e$message))
  })
  
  if (startsWith(status[[1]], "Error:")) {
    abort("Cas-OFFinder not installed properly! Please refer to https://github.com/snugel/cas-offinder for installation.")
  } else {
    # inform("Cas-OFFinder installed properly!")
    return(ifelse(any(grepl("Type: GPU", status)), "G", "C"))
  }
}

# Prepare bulge wrapper
prepCasOFFinderBulge <- function(version = c("2.4.1", "3.0.0b3")) {
  if (version == "2.4.1") {
    wrapper_path <- system.file("extdata", "Cas-OFFinder/cas-offinder-bulge", package = "CRISPRseek")
    
    # Check and install Python if needed
    py_config()
    if (py_available()) {
      pyv <- unlist(strsplit(py_discover_config()$version, ".", fixed = TRUE))
      if (length(pyv) < 2 || !( (pyv[1] == 3 && pyv[2] >= 5) )) {
        abort("Please install Python 3.5 or higher for Cas-OFFinder!")
      }
    } else {
      tryCatch(
        install_miniconda(),
        error = function(e) { print(e) }
      )
      py_config()
    }
    
    # Load cas_offinder_bulge()
    source_python(wrapper_path)
    return(cas_offinder_bulge)
  } else if (version == "3.0.0b3") {
    cas_offinder <- function(input_txt, device, output_name, binary_path) {
      cmd <- paste(paste0('"', binary_path, '"'), paste0('"', input_txt, '"'), device, paste0('"', output_name, '"'), "> /dev/null 2>&1", sep = " ")
      tryCatch({
        system(cmd)
      }, error = function(e) {
        inform(paste0("Error:", e$message))
      })
    }
    return(cas_offinder)
  } else {
    abort("Unsupported CasOFFinder version!")
  }
}

# Prepare CasOFFinder input files: genome
prepCasOFFinderInputGenome <- function(BSgenomeName, genomeSeqFile, chromToSearch = "all", chromToExclude = NULL, outdir) {
  if (!is.null(genomeSeqFile) && file.exists(genomeSeqFile)) {
    ref <- readDNAStringSet(genomeSeqFile)
    chrs <- names(ref)
  } else {
    ref <- BSgenomeName
    chrs <- seqnames(ref)
  }
  outdir <- tempfile(pattern = "ref_", tmpdir = outdir, fileext = "")
  dir.create(outdir)
  lapply(chrs, function(chr) {
    if (chromToSearch[[1]] == "all" || chr %in% chromToSearch && !(chr %in% chromToExclude)) {
      seq <- DNAStringSet(ref[[chr]])
      names(seq) <- chr
      fa <- file.path(outdir, paste0(chr, ".fasta"))
      writeXStringSet(seq, filepath = fa, format = "fasta")
    }
  }
)
  return(outdir)
}

# Prepare CasOFFinder input files: input.txt 
prepCasOFFinderInputTxt <- function(gRNA_PAM, 
                                    PAM.size = 3, 
                                    PAM.pattern = "NNG$|NGN$", 
                                    PAM.location = "3prime", 
                                    max.mismatch = 3, 
                                    DNA_bulge = 2, 
                                    RNA_bulge = 2,
                                    chr_path,
                                    cas_offinder_version,
                                    outdir = tempdir()) {
  # input.txt format: line1 is chr fa dir path; line2 is gRNA(all N)+PAM.pattern, DNA bulges, RNA bulges; line3 and beyond are gRNA + PAM(N) and allowed mismatches 
  # /Users/kaihu/Downloads/chr
  # NNNNNNNNNNNNNNNNNNNNNRG 3 3
  # GGCCGACCTGTCGCTGACGCNNN 3
  
  line1 <- chr_path 
  
  line2_backbone <- paste0(rep("N", gRNA_PAM[[1]]@length - PAM.size), collapse = "")
  pattern <- translatePAMpattern(PAM.pattern, direction = "toIUPAC")
  v_pattern <- unlist(strsplit(pattern, "\\|"))
  line2s <- sapply(v_pattern, function(pattern) {
    if (grepl("\\$", pattern)) {
      line2 <- paste0(line2_backbone, gsub("\\$", "", pattern))
    } else if (grepl("\\^", pattern)) {
      line2 <- paste0(gsub("\\^", "", pattern), line2_backbone)
    } else { abort("Neither $ or ^ exists in PAM.pattern!") }
    paste(line2, DNA_bulge, RNA_bulge, collapse = " ")
  })
  
  line_gRNAs <- sapply(names(gRNA_PAM), function(g_name) {
    g <- gRNA_PAM[[g_name]]
    if (PAM.location == "3prime") {
      line <- paste0(g[1:(length(g) - PAM.size)], paste0(rep("N", PAM.size), collapse = ""))
    } else if (PAM.location == "5prime") {
      line <- paste0(paste0(rep("N", PAM.size), collapse = ""), g[(PAM.size+1):length(g)])
    } else { abort("Wrong PAM.location!") }
    
    if (cas_offinder_version == "2.4.1") {
      paste(line, max.mismatch, collapse = " ")
    } else if (cas_offinder_version == "3.0.0b3") {
      paste(line, max.mismatch, g_name, collapse = " ")
    }
  })
  
  input_txt <- lapply(line2s, function(line2) {
    unname(c(line1, line2, line_gRNAs))
  })
  
  input_names <- file.path(outdir, paste0("input", seq_along(input_txt), "_", basename(tempfile(pattern = "", fileext = ".txt"))))
  
  invisible(mapply(function(v, v_name) {
    writeLines(v, con = v_name)
  }, input_txt, input_names))
  
  return(input_names)
}

# Filter results obtained from CasOFFinder v3.0.0b3
filterCasOffinder <- function(df_bulge = NULL,
                              max.mismatch = NULL, 
                              gRNA.size = NULL,
                              PAM = NULL,
                              PAM.size = NULL,
                              PAM.pattern = NULL, 
                              PAM.location = NULL,
                              allowed.mismatch.PAM = NULL, 
                              baseEditing = NULL,
                              targetBase = NULL,
                              editingWindow = NULL) {
  # 1. Ensure the max.mismatch to offtarget DNA is correct (only considering the gRNA part excluding PAM per CRISPRseek definition)
  # The results from Cas-OFFinder is clearly not right regarding the max.mismatch column
  df_bulge$mismatches_c <- str_count(df_bulge$DNA, "[a-z]")
  df_bulge <- df_bulge[df_bulge$mismatches_c <= max.mismatch, ]
  
  # 2. Ensure offtarget DNA contains the PAM.pattern: Cas-OFFinder results should all contain PAM.pattern, just a safety layer in case
  PAM.pattern <- translatePAMpattern(PAM.pattern, direction = "fromIUPAC")
  df_bulge <- df_bulge[grepl(PAM.pattern, df_bulge$DNA, perl = TRUE), ]
  
  # 3. Ensure allowed.mismatch.PAM: by comparing offtarget DNA PAM with PAM (canonical PAM)
  if (PAM.location == "3prime") {
    PAM_sequence <- substr(df_bulge$DNA, nchar(df_bulge$DNA) - PAM.size + 1, nchar(df_bulge$DNA))
  } else if (PAM.location == "5prime") {
    PAM_sequence <- substr(df_bulge$DNA, 1, PAM.size)
  }
  
  canonical_PAM <- DNAString(PAM)
  n_PAM_mismatch <- unlist(lapply(DNAStringSet(PAM_sequence), function(PAM_s) {
    neditAt(PAM_s, canonical_PAM, fixed = FALSE)
  }))
  
  df_bulge <- df_bulge[n_PAM_mismatch <= allowed.mismatch.PAM, ]
  
  # 4. Deal with Base editing requirement
  if (baseEditing) {
    editingSeq <- str_sub(as.character(DNAStringSet(df_bulge$DNA)), min(editingWindow), max(editingWindow))
    num_target_base <- str_count(editingSeq, toupper(targetBase))
    df_bulge <- df_bulge[num_target_base > 0, ]
  }
  
  return(df_bulge)
}

# Format results obtained from CasOFFinder after filtering    
# strand, chrom, chromStart, chromEnd, name, gRNAplusPAM, gRNAPlusPAM_bulge, OffTargetSequence_bulge, n.mismatch, n.DNABulge, n.RNABulge, forViewInUCSC, score
formatCasOffinder <- function(df_bulge = NULL,
                              gRNAs = NULL) {
  chrom = df_bulge$chr
  chromStart = df_bulge$start_0_based + 1
  chromEnd = df_bulge$start_0_based + nchar(gsub("-", "", df_bulge$DNA))
  hits_bulge <- data.frame(strand = df_bulge$strand,
                           chrom = chrom, 
                           chromStart = chromStart, 
                           chromEnd = chromEnd, 
                           name = df_bulge$gRNA_id,  
                           gRNAPlusPAM = as.character(gRNAs[df_bulge$gRNA_id]),
                           gRNAPlusPAM_bulge = df_bulge$gRNA,
                           OffTargetSequence_bulge = df_bulge$DNA,
                           n.mismatch = df_bulge$mismatches_c,
                           n.RNABulge = ifelse(df_bulge$bulge_type == "RNA", df_bulge$bulge_size, 0),
                           n.DNABulge = ifelse(df_bulge$bulge_type == "DNA", df_bulge$bulge_size, 0),
                           forViewInUCSC = paste0(chrom, ":", chromStart, "-", chromEnd),
                           score = rep(100, nrow(df_bulge)))
  return(hits_bulge)
}


# Add columns required by buildFeatureVectorForScoringBulge to hits
addAlnInfoToHits <- function(hits = NULL,
                             gRNA.size = NULL,
                             PAM.size = NULL,
                             PAM.location = NULL) {
  # add offTarget_sequence
  hits$offTarget_sequence <- str_to_upper(str_replace_all(hits$OffTargetSequence_bulge, "-", ""))
  
  # determine PAM.sequence
  if (PAM.location == "3prime") {
    hits$PAM.sequence <- str_sub(hits$OffTargetSequence_bulge, -PAM.size)
  } else if (PAM.location == "5prime") {
    hits$PAM.sequence <- str_sub(hits$OffTargetSequence_bulge, 1, PAM.size)
  }
  
  # determine pos.mismatch: in gRNA coordinates, where mismatches occur
  replace_uppercase <- function(str1, str2, pattern = "[a-z]") {
    ## given str2, find pos of pattern, and replace same pos in str1 (str1 all have uppercase) into lower case letters
    lower_positions <- which(str_detect(str_split(str2, "")[[1]], pattern))
    chars1 <- str_split(str1, "")[[1]]
    chars1[lower_positions] <- str_to_lower(chars1[lower_positions])
    return(str_c(chars1, collapse = ""))
  }
  gRNAPlusPAM_bulge_mismatch <- mapply(replace_uppercase, 
                                       hits$gRNAPlusPAM_bulge, 
                                       hits$OffTargetSequence_bulge, 
                                       pattern = "[a-z]")
  hits$pos.mismatch <- lapply(gRNAPlusPAM_bulge_mismatch, function(str) {
    str <- str_replace_all(str, "-", "")
    which(str_detect(str_split(str, "")[[1]], "[a-z]"))
  })
  
  # determine pos.insertion: in gRNA coordinates, where offtargets have insertion (DNA bulge)
  ## when determining the pos.insertion, use gRNA coordinates: count number of non-hyphen nucleotide before each hyphen, for consecutive hyphens, use continous numbering but don't allow larger than gRNA.size positions. Caveat: the CFD algorithm may not support multiple bulges.
  count_non_hypen_letters <- function(str, pattern = "-") {
    count <- 0
    result <- c()
    for (char in strsplit(str, "")[[1]]) {
      if (char == pattern) {
        if (length(result) > 0) {
          if (count == result[[length(result)]]) {
            if (count <  gRNA.size) {
              result <- c(result, count + 1)
            }
          } else {
            result <- c(result, count)
          }
        } else {
          result <- c(result, count)
        }
      } else if (char != "" && char != "-") {
        count <- count + 1
      }
    }
    return(result)
  }
  hits$pos.insertion <- lapply(hits$gRNAPlusPAM_bulge, count_non_hypen_letters, pattern = "-")
  
  # determine pos.deletion: in gRNA coordinates, where offtargets have deletion (RNA bulge)
  gRNAPlusPAM_bulge_DNA_bulge <- mapply(replace_uppercase, 
                                        hits$gRNAPlusPAM_bulge, 
                                        hits$OffTargetSequence_bulge, 
                                        pattern = "-")
  hits$pos.deletion <- lapply(gRNAPlusPAM_bulge_DNA_bulge, function(str) {
    str <- str_replace_all(str, "-", "")
    which(str_detect(str_split(str, "")[[1]], "[a-z]"))
  })
  
  # determine guideAlignment2OffTarget: gRNA alignment to offtarget seq
  getAlnToOffTarget <- function(gRNA = NULL, 
                                offtarget = NULL, 
                                PAM.location = NULL,
                                PAM.size = NULL) {
    if (PAM.location == "3prime") {
      gRNA <- str_sub(gRNA, 1, -PAM.size-1)
      offtarget <- str_sub(offtarget, 1, -PAM.size - 1)
    } else if (PAM.location == "5prime") {
      gRNA <- str_sub(gRNA, PAM.size + 1)
      offtarget <- str_sub(offtarget, PAM.size + 1)
    }
    
    # initialize alignment with all ., meaning matches
    res <- rep(".", nchar(gRNA))
    
    # if mismatch, use offtarget seq base:
    gRNA_bases <- str_split(str_to_upper(gRNA), "")[[1]]
    offtarget_bases <- str_split(str_to_upper(offtarget), "")[[1]]
    mismatch_pos <- which(gRNA_bases != offtarget_bases)
    mismatch_pos <- mismatch_pos[offtarget_bases[mismatch_pos] != "-"]
    res[mismatch_pos] <- offtarget_bases[mismatch_pos]
    
    # replace insertions in DNA (DNA bulge) with "^":
    dna_bulge_pos <- which(str_detect(str_split(gRNA, "")[[1]], "-"))
    res[dna_bulge_pos] <- "^"
    
    # replace deletions in DNA (RNA bulge) with "-":
    rna_bulge_pos <- which(str_detect(str_split(offtarget, "")[[1]], "-"))
    res[rna_bulge_pos] <- "-"
    
    guideAlignment2OffTarget <- str_c(res, collapse = "")
    return(guideAlignment2OffTarget)
  }
  
  hits$guideAlignment2OffTarget <- mapply(getAlnToOffTarget, 
                                          hits$gRNAPlusPAM_bulge, 
                                          hits$OffTargetSequence_bulge, 
                                          PAM.location = PAM.location, 
                                          PAM.size = PAM.size)
  
  return(hits)
}


# Get off-target score with bulges:
.updateScore <- function(fv.geThan1, 
                         fv.lessThan1,
                         col.prefix = "IsInsertion.pos",
                         weights, 
                         position, 
                         type) {
  if (dim(fv.geThan1)[1] > 0) {
    indel.pos <- fv.geThan1[,grep(col.prefix, colnames(fv.geThan1))]
    indel.pos <- apply(indel.pos, 1, as.numeric)
    pos <- grep(col.prefix, colnames(fv.geThan1))
    min.pos <- min(pos)

    if (length(pos) > 0) {
      score.new <- unlist(lapply(1:nrow(fv.geThan1), function(i) {
        indel.index <- pos[fv.geThan1[i, pos] == 1] - min.pos + 1
        if (col.prefix == "IsInsertion.pos") {
          this.type <- unlist(strsplit(as.character(fv.geThan1[i,]$gRNA.insertion), ","))
        } else {
          this.type <- unlist(strsplit(as.character(fv.geThan1[i,]$gRNA.deletion), ","))
        }
        
        score.new1 <- fv.geThan1[i, ]$score
        for (j in 1:length(indel.index)) {
          score.new1 <- score.new1 * weights[position == indel.index[j] & type == this.type[j]]
        }
        
        return(score.new1)
      }))
      fv.geThan1$score <- score.new
    }
  }
  
  if (dim(fv.geThan1)[1] > 0) {
    score <- fv.geThan1
    if (dim(fv.lessThan1)[1] > 0) {
      score <- rbind(fv.lessThan1, fv.geThan1)
    }
  } else {
    score <- fv.lessThan1
  }
  
  return(score)
}

#' @noRd
#' @title Obtain scores for off-targets with bulges
#' @importFrom BiocGenerics subset unlist lapply rbind
#' @importFrom hash hash values
#' @importFrom rio import
#' @importFrom hash hash
#' @examples
#' 
#' getOfftargetScoreBulge(featureVectors, 
#'                        mismatch.activity.file = mismatch.activity.file_default_xlsx())
getOfftargetScoreBulge <- function(featureVectors,
                                   subPAM.activity = subPAM.activity_default,
                                   mismatch.activity.file = mismatch.activity.file_default_xlsx()) {
  featureVectors <- getOfftargetScore2(featureVectors)
  deletion.activity <- import(mismatch.activity.file, which = 2)
  insertion.activity <- import(mismatch.activity.file, which = 3)
  required.col.ins <- c("Insertion.Type", "Position", "Percent.Active")
  required.col.del <- c("Deletion.Type", "Position", "Percent.Active")
  
  if (length(intersect(colnames(insertion.activity), required.col.ins)) != length(required.col.ins)) {
    stop("Please rename the insertion sheet of the activity file column to contain at least these 3 column names: Insertion.Type, Position, Percent.Active\n")
  }
  
  if (length(intersect(colnames(deletion.activity), required.col.del)) != length(required.col.del)) {
    stop("Please rename the deletion sheet of the activity file column to contain at least these 3 column names: Deletion.Type, Position, Percent.Active\n")
  }
  
  colnames(insertion.activity)[colnames(insertion.activity) == "Insertion.Type"] <- "Type"
  colnames(deletion.activity)[colnames(deletion.activity) == "Deletion.Type"] <- "Type"
  
  #mismatch.activity[mismatch.activity$Mismatch.Type == "rA:dG" &
  #    mismatch.activity$Position == 10,]$Percent.Active
  ##### by default weights is a column vector
  ##### the mismatch activity  is given as pos 20, 19, 18,....1 distance from PAM,    
  ##### Position named as 1, 2, 3, .....20 though
  ##### and the featureVectors is in the same order now
  ##### so no need to reverse any more. weights = rev(weights)
  
  # Infer n.insertion (equals length of pos.insertion)
  # Infer n.deletion (equals to the length of pos.deletion)
  featureVectors$n.insertion <- lapply(featureVectors$pos.insertion, function(x) length(x))
  featureVectors$n.deletion <- lapply(featureVectors$pos., function(x) length(x))
  
  fv.geThan1 <- subset(featureVectors, as.numeric(as.character(featureVectors$n.insertion)) >= 1)
  fv.lessThan1 <- subset(featureVectors, as.numeric(as.character(featureVectors$n.insertion)) < 1)
  
  featureVectors <- .updateScore(fv.geThan1, 
                                 fv.lessThan1,
                                 col.prefix = "IsInsertion.pos",
                                 weights = insertion.activity$Percent.Active,
                                 type = insertion.activity$Type,
                                 position = insertion.activity$Position)
  
  # do the same for deletion
  fv.geThan1 <- subset(featureVectors, as.numeric(as.character(featureVectors$n.deletion)) >= 1)
  fv.lessThan1 <- subset(featureVectors, as.numeric(as.character(featureVectors$n.deletion)) < 1)
  
  score <- .updateScore(fv.geThan1, 
                        fv.lessThan1,
                        col.prefix = "IsDeletion.pos",
                        weights = deletion.activity$Percent.Active,
                        type = deletion.activity$Type,
                        position = deletion.activity$Position)
  
  score$alignment <- as.character(score$alignment)
  score$score <- round(score$score, 6)
  score  <- score[order(c(score$name,score$score),decreasing = TRUE), ]
  
  return(unique(score[!is.na(score$score), ]))
}



