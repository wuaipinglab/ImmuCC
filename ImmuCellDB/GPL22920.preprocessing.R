#!/usr/bin/Rscript --slave
#### Description:
     # Convert Raw microarray files profiled in Human Genome Hglue2.0 platform into expression values

################################################################################################################################################
###                                         Data Processing for Hglue2.0
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
  library(hglue20hsentrezgcdf)
  library(hglue20hsentrezg.db)

  Data <- ReadAffy(cdfname = "hglue20hsentrezgcdf")
  cat("Raw CEL has been successfully read into AffyBatch object!\n")
  saveRDS(Data, file=paste(Result.dir, Result.name, ".AffyData.customCDF.RDS", sep=""), compress=F)

  # rma normalization!
  eset <- rma(Data)
  cat("RMA normalization has accompanished!\n")
  saveRDS(eset, file=paste(Result.dir, Result.name, ".esetRMA.customCDF.RDS", sep=""), compress=F)

  ###########################################################
  ID <- featureNames(eset)
  GS <- as.matrix(getSYMBOL(ID, 'hglue20hsentrezg.db'))    # Mapping the Probeset into the correspondent Gene Symbol
  EG <- as.matrix(getEG(ID, 'hglue20hsentrezg.db'))        # Mapping the Probeset into the correspondent EntrezGene

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

