#Sys.setenv(PATH = paste("/Users/ZHUJ/anaconda2/bin", Sys.getenv("PATH"), sep=":"))
#system("python --version")
#Python 2.7.15 :: Anaconda, Inc.

#extendedSequence="TCTGGTGTCCTCCACACCAGAATCAGGGGT"

#CRISPRseek:::calculategRNAEfficiency2(extendedSequence,aa.cut = -1, per.peptide = -1)

#' @importFrom reticulate py_discover_config
#' @importFrom reticulate conda_list
#' @importFrom reticulate conda_create
#' @importFrom reticulate conda_install
#' @importFrom reticulate py_run_string
#' @importFrom BiocGenerics lapply

calculategRNAEfficiency2 <- function(extendedSequence, aa.cut = -1, per.peptide = -1) {
    system2("python", args = "--version", stderr="pythonVersion.txt")
    #if(grep( "2.7", read.table("pythonVersion.txt", sep="", header=FALSE)[1,2]))
    #The grep("2.7") method above is not robust because version can ben NULL and grep can't handle NULL properly.

    ### Dev notes start ###
    ## First check if there is a system Python2.7 that is equipped with all required modules, if so, use that, no installation is required.
    ## Otherwise, need to create py2 env with reticulate:
    ##    First check if py2 exists and all modules are installed: if not, recreate py2.
    ### Dev notes end ###

    # Check if system Python2.7 pyenv is satisfied or not:
    system_pyenv = FALSE # indicate if system_pyenv is ready or not
    # First, delete existing system_py2_missing_modules.txt file:
    if (file.exists("system_py2_missing_modules.txt")) {
      file.remove("system_py2_missing_modules.txt")
    }
    if (py_discover_config()$version == "2.7") {
      py_run_string("import sys")
      py_run_string("import subprocess")
      py_run_string("import pkg_resources")
      py_run_string("required = {'numpy', 'scipy', 'scikit-learn', 'pandas', 'matplotlib', 'biopython'}")
      py_run_string("installed = {pkg.key for pkg in pkg_resources.working_set}")
      py_run_string("missing = required - installed")
      py_run_string("if len(missing) > 0: \n\twith open('system_py2_missing_modules.txt', 'w') as fp: \n\t\tfp.write(str(missing))")

      if (!(file.exists("system_py2_missing_modules.txt"))) {
        system_pyenv = TRUE
      }
    }

    # Setup conda pyenv if not system Python not satisfied:
    if (!(system_pyenv)) {
      cat("System python2.7 env not ready, preparing conda (py2) ...\n")
      tryCatch(
        expr = {
          # Check if py2 exists:
          # First, delete existing py2_missing_modules.txt file:
          if (file.exists("py2_missing_modules.txt")) {
            file.remove("py2_missing_modules.txt")
          }
          path_py2 <- conda_list()[which(conda_list()$name == "py2")[1], 2]
          bin_py2 <- substring(path_py2, 1, nchar(path_py2) - 7)
          py2 <- FALSE

          if (!(nchar(path_py2) == 0)) {
            cat("CHECK: py2 env already created. Checking required modules ...\n")

            Sys.setenv(RETICULATE_PYTHON = path_py2)
            Sys.setenv(PATH = paste(bin_py2, Sys.getenv("PATH"), sep = ":"))
            py_run_string("import sys")
            py_run_string("import subprocess")
            py_run_string("import pkg_resources")
            py_run_string("required = {'numpy', 'scipy', 'scikit-learn', 'pandas', 'matplotlib', 'biopython'}")
            py_run_string("installed = {pkg.key for pkg in pkg_resources.working_set}")
            py_run_string("missing = required - installed")
            py_run_string("if len(missing) > 0: \n\twith open('py2_missing_modules.txt', 'w') as fp: \n\t\tfp.write(str(missing))")

            if (!(file.exists("py2_missing_modules.txt"))) {
              py2 <- TRUE
            } else {
              cat("CHECK: py2 env missing some required modules, need to recreate ...\n")
            }
          } else {
            cat("CHECK: py2 env not created yet. Need to create ...\n")
          }

          if (py2) {
            cat("CHECK: py2 env and all required modules all satisfied.\n")
          } else {
            cat("CHECK: creating py2 ...\n")
            conda_create("py2", python="2.7.15")
            conda_install("py2", "numpy")
            conda_install("py2", "pandas")
            conda_install("py2", "scipy", pip = TRUE)
            conda_install("py2", "scikit-learn==0.16.1", pip = TRUE)
            conda_install("py2", "matplotlib", pip = TRUE)
            conda_install("py2", "biopython") # pip will install 1.77, which is not supported by python2.7 any longer, conda will install 1.76
            path_py2 <- conda_list()[which(conda_list()$name == "py2")[1], 2] # first here in case more than two py2 avail, conda_install() seems would always install to the first one.
            bin_py2 <- substring(path_py2, 1, nchar(path_py2) - 7)
            Sys.setenv(RETICULATE_PYTHON = path_py2)
            Sys.setenv(PATH = paste(bin_py2, Sys.getenv("PATH"), sep = ":"))
            cat("Rule_Set_2 python2.7 env (py2) created!\n")
          }
        },
        error = function(e) {
          print(e)
          stop("Python 2.7 is required for calculating gRNA efficacy using Rule_Set_2_scoring published in JG. Doench, et al., Nbt, Jan 2016! Check CRISPRseek manual for details.")
        }
      )
    } else {
      cat("Rule_Set_2 python2.7 env satisfied, continuing ...\n")
    }

    # Run rule_set_2:
    origDir <- getwd()
    pythonDir <- system.file("extdata/Rule_Set_2_scoring_v1/analysis/",package = "CRISPRseek")
    if (dir.exists(pythonDir)) {
      setwd(pythonDir)
      tryCatch(
        (efficiency <- unlist(lapply(extendedSequence, function(thisSeq) {
          py_command <- paste("python rs2_score_calculator.py -seq ", thisSeq, " -aa-cut ", aa.cut, " -per-peptide ", per.peptide, sep="")
          system(py_command, intern=TRUE)
        }))),
        error = function(e) {
          print(e);
        }
      )
      efficiency <- unlist(lapply(efficiency, function(temp) {strsplit(temp, ": ")[[1]][2]}))
      setwd(origDir)
    } else {
      stop("Rule_Set_2_scoring_v1 python code not found in the extdata directory of CRISPRseek, please update CRISPRseek to the newest version first!")
    }
   efficiency
}
