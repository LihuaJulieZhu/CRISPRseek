foldgRNAs <- function(gRNAs.withoutPAM, 
    gRNA.backbone="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
    temperature = 37)
{
    GeneRfold.installed <- 1 
    tryCatch({library(GeneRfold)}, error = function(e) 
     {
	print(paste("MY_ERROR:  ",e))	
	GeneRfold.installed <<- 0	
     })
     if ( GeneRfold.installed ) { 
       fd1 <- GeneRfold::fold(s = gRNA.backbone, t = temperature)	
       mfe.backbone <- round(fd1[[2]], 2)
       gRNAs <- gsub("T", "U", gRNAs.withoutPAM)
       fd2 <- matrix(nrow = length(gRNAs), ncol = 2)
       #fd2 <- do.call(rbind, lapply( gRNAs, function(s) {
       #   unlist(fold(paste(s, gRNA.backbone, sep=""), t=temperature))}))
       for ( i in 1:length(gRNAs))
       {
          fd2[i,] <- unlist(GeneRfold::fold(paste(gRNAs[i], gRNA.backbone, sep=""), 
             t=temperature)) 
       }
       colnames(fd2) <- c("bracket.notation", "mfe.sgRNA")
       fd2[,2] <- round(as.numeric(unlist(fd2[,2])))
       mfe.diff <- round(as.numeric(unlist(fd2[,2])) - mfe.backbone,2)
       if (length(gRNAs.withoutPAM) > 1)
       	   fd <- cbind(gRNAs.withoutPAM, fd2, mfe.diff,  mfe.backbone)
       else
	{
           fd <- t(data.frame(c(gRNAs.withoutPAM, fd2, mfe.diff, mfe.backbone)))
	   colnames(fd) <- c("gRNAs.withoutPAM", "bracket.notation", "mfe.sgRNA",
               "mfe.diff",  "mfe.backbone")
	}
       fd
   }
   else
   { 
       fd <- gRNAs.withoutPAM 
   }
}
