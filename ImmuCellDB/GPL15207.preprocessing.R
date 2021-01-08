
#!/usr/bin/Rscript --slave
#### Description:
     # Convert Raw microarray files profiled in Human Genome 133platforms into expression values!

# /gluster/home/chenziyi/software/R/R-3.2.5/bin/R

################################################################################################################################################
###                                         Data Processing for Hgu133a2.0
################################################################################################################################################

  cat("********************************************************************************************************************************************\n")
  argv <- commandArgs(TRUE)
  options(stringsAsFactors=F)

  Data.dir <- argv[1]
  Result.name <- argv[2]
  Result.dir <- argv[3]
  setwd(Data.dir)
  cat(Result.name, "/n")

  library(affy)
  library(annotate)
  library(primeviewhsentrezgcdf)
  library(primeviewhsentrezg.db)

  # Performs the Wilcoxon signed rank-based gene expression presence/absence detection algorithm!
  # first implemented in the Affymetrix Microarray Suite version 5,The mas5calls method for AffyBatch \n
  # returns an ExpressionSet with  calls accessible with exprs(obj) and p-values available with assayData(obj)
  PMA_calls <- function(affydata){
      # Argument:
      #     affydata:raw file produced by ReadAffy()
      #     PMA_file: name of result
      cat("Starting to calculate PMA calls using mas5calls function......", "\n")       
      calls <- mas5calls(affydata) 
      calls <- exprs(calls)
      stat <- apply(calls, 1, function(x){calls.stat <- table(x); calls.prop <- calls.stat/sum(calls.stat); calls.prop})
      present <- unlist(lapply(stat, function(x){x["P"]}))
      marginal <- unlist(lapply(stat, function(x){x["M"]}))
      absent <- unlist(lapply(stat, function(x){x["A"]}))
      pma <- data.frame(Present=present,
                        Absent=absent,
                        Marginal=marginal
                       )
      pma[is.na(pma)==T] <- 0
      rownames(pma) <- names(stat)
      pma
  }

  Data <- ReadAffy(cdfname = "primeviewhsentrezgcdf")
  cat("Raw CEL has been successfully read into AffyBatch object!\n")
  saveRDS(Data, file=paste(Result.dir, Result.name, ".AffyData.customCDF.RDS", sep=""), compress=F)

  # rma normalization!
  eset <- rma(Data)
  cat("RMA normalization has accompanished!\n")
  saveRDS(eset, file=paste(Result.dir, Result.name, ".esetRMA.customCDF.RDS", sep=""), compress=F)

  ###########################################################
  ID <- featureNames(eset)
  GS <- as.matrix(getSYMBOL(ID, 'primeviewhsentrezg.db'))    # Mapping the Probeset into the correspondent Gene Symbol
  EG <- as.matrix(getEG(ID, 'primeviewhsentrezg.db'))        # Mapping the Probeset into the correspondent EntrezGene

  expression <- exprs(eset)
  colnames(expression)<- sub("_.*", "", colnames(expression))
  colnames(expression) <- sub("\\..*", "", colnames(expression))
  cat("Expression data stored in eSets objects has been successfully accessed!\n")

  cols <-  c("ProbeSet", "EntrezGene", "GeneSymbol")
  mapping <- cbind(rownames(expression), EG, GS)
  mapping <- mapping[which(mapping[, 3] != "NA"), ]            # remove NAs
  mapping <- mapping[order(mapping[, 1]), ]                    # sort by gene name 
  colnames(mapping) <- cols
  rownames(mapping) <- mapping[, 1]

  expression.knownGene <- expression[rownames(mapping), ]                      # expression matrix with Gene Symbols
  rownames(expression.knownGene) <- mapping[rownames(expression.knownGene), 3]

  saveRDS(mapping, expression, expression.knownGene, file=paste(Result.dir, Result.name, ".FinalExpressionRMA.RDS", sep=""), compress=F)
  write.table(expression.knownGene, file=paste(Result.dir, Result.name, ".ExpressionArray.customCDF.txt", sep=""), sep="\t", col.names=F, row.names=F, quote=FALSE)

