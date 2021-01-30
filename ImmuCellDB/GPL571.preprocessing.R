
#!/usr/bin/Rscript --slave
#### Description:
     # Convert Raw microarray files profiled in Huu133a2.0 platforms into expression values!

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
  library(frma)
  library(annotate)
  library(hgu133a2hsentrezgcdf)
  library(hgu133a2hsentrezg.db)
  library(hgu133a2frmavecs)

  Data <- ReadAffy(cdfname = "hgu133a2hsentrezgcdf")
  cat("Raw CEL has been successfully read into AffyBatch object!\n")
  saveRDS(Data, file=paste(Result.dir, Result.name, ".AffyData.customCDF.RDS", sep=""), compress=F)

  # Frma normalization!
  eset <- frma(Data, normalize="quantile", output.param = "hgu133a2hsentrezgfrmavecs")
  cat("FRMA normalization has accompanished!\n")
  saveRDS(eset, file=paste(Result.dir, Result.name, ".esetFRMA.customCDF.RDS", sep=""), compress=F)

  ###########################################################

  # Calculation the global normalized unscaled standard error(GNUSE)
  gnuse <- GNUSE(eset, type = "stat")
  colnames(gnuse)<- gsub("_.*", "", colnames(gnuse))
  colnames(gnuse) <- gsub("\\..*", "", colnames(gnuse))

  gnusecutoff <- gnuse[1, ] <= 1.25
  gnuse.selected <- gnuse[, gnusecutoff]
  colnames(gnuse.selected)<- gsub("_.*", "", colnames(gnuse.selected))
  colnames(gnuse.selected) <- gsub("\\..*", "", colnames(gnuse.selected))
  cat("GNUSE has calculated!\n")
  write.csv(t(gnuse.selected), file=paste(Result.dir, Result.name, ".gnuse_quatified.csv", sep=""))  
  saveRDS(gnuse, gnuse.selected, file = paste(Result.dir, Result.name, ".GNUSE.RDS", sep=""), compress=F)

  ID <- featureNames(eset)
  GS <- as.matrix(getSYMBOL(ID, 'hgu133a2hsentrezg.db'))    # Mapping the Probeset into the correspondent Gene Symbol
  EG <- as.matrix(getEG(ID, 'hgu133a2hsentrezg.db'))        # Mapping the Probeset into the correspondent EntrezGene

  expression <- exprs(eset)
  colnames(expression)<- gsub("_.*", "", colnames(expression))
  colnames(expression) <- gsub("\\..*", "", colnames(expression))
  expression <- expression[, colnames(gnuse.selected)]
  cat("Expression data stored in eSets objects has been successfully accessed!\n")

  cols <-  c("ProbeSet", "EntrezGene", "GeneSymbol")
  mapping <- cbind(rownames(expression), EG, GS)
  mapping <- mapping[which(mapping[, 3] != "NA"), ]            # remove NAs
  mapping <- mapping[order(mapping[, 1]), ]                    # sort by gene name 
  colnames(mapping) <- cols
  rownames(mapping) <- mapping[, 1]

  expression.knownGene <- expression[rownames(mapping), ]                      # expression matrix with Gene Symbols
  rownames(expression.knownGene) <- mapping[rownames(expression.knownGene), 3]

  saveRDS(mapping, expression, expression.knownGene, file=paste(Result.dir, Result.name, ".FinalExpressionFRMA.RDS", sep=""), compress=F)
  write.table(expression.knownGene, file=paste(Result.dir, Result.name, ".ExpressionArray.customCDF.txt", sep=""), sep="\t", col.names=F, row.names=F, quote=FALSE)

