
#!/usr/bin/Rscript --slave
#### Description:
     # Convert Raw microarray files profiled in Human Genome Hthgu133a platform into expression values!

################################################################################################################################################
###                                         Data Processing for Ht hgu133a
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
  library(hthgu133ahsentrezgcdf)
  library(hthgu133ahsentrezg.db)
  library(hthgu133afrmavecs)

  Data <- ReadAffy(cdfname = "hthgu133ahsentrezgcdf")
  cat("Raw CEL has been successfully read into AffyBatch object!\n")
  save(Data, file=paste(Result.dir, Result.name, ".AffyData.customCDF.RData", sep=""))

  # Frma normalization!
  eset <- frma(Data, normalize="quantile", output.param = "hthgu133ahsentrezgfrmavecs")
  cat("FRMA normalization has accompanished!\n")
  save(eset, file=paste(Result.dir, Result.name, ".esetFRMA.customCDF.RData", sep=""))

  ###########################################################

  # Calculation the global normalized unscaled standard error(GNUSE)
  gnuse <- GNUSE(eset, type = "stat")
  colnames(gnuse)<- sub("_.*", "", colnames(gnuse))
  colnames(gnuse) <- sub("\\..*", "", colnames(gnuse))

  gnusecutoff <- gnuse[1, ] <= 1.25
  gnuse.selected <- gnuse[, gnusecutoff]
  colnames(gnuse.selected)<- gsub("_.*", "", colnames(gnuse.selected))
  colnames(gnuse.selected) <- gsub("\\..*", "", colnames(gnuse.selected))
  cat("GNUSE has calculated!\n")
  write.csv(t(gnuse.selected), file=paste(Result.dir, Result.name, ".gnuse_quatified.csv", sep=""))  
  save(gnuse, gnuse.selected, file = paste(Result.dir, Result.name, ".GNUSE.RData", sep=""))

  ID <- featureNames(eset)
  GS <- as.matrix(getSYMBOL(ID, 'hthgu133ahsentrezg.db'))    # Mapping the Probeset into the correspondent Gene Symbol
  EG <- as.matrix(getEG(ID, 'hthgu133ahsentrezg.db'))        # Mapping the Probeset into the correspondent EntrezGene

  expression <- exprs(eset)
  colnames(expression)<- sub("_.*", "", colnames(expression))
  colnames(expression) <- sub("\\..*", "", colnames(expression))
  expression <- expression[, colnames(gnuse.selected)]
  cat("Expression data stored in eSets objects has been successfully accessed!\n")
  save(expression, file=paste(Result.dir, Result.name, ".Expression.RData", sep=""))

  cols <-  c("ProbeSet", "EntrezGene", "GeneSymbol")
  mapping <- cbind(rownames(expression), EG, GS)
  mapping <- mapping[which(mapping[, 3] != "NA"), ]            # remove NAs
  mapping <- mapping[order(mapping[, 1]), ]                    # sort by gene name 
  colnames(mapping) <- cols
  rownames(mapping) <- mapping[, 1]

  expression.knownGene <- expression[rownames(mapping), ]                      # expression matrix with Gene Symbols
  rownames(expression.knownGene) <- mapping[rownames(expression.knownGene), 3]

  save(mapping, expression, expression.knownGene, file=paste(Result.dir, Result.name, ".FinalExpressionFRMA.RData", sep=""))
  write.table(expression.knownGene, file=paste(Result.dir, Result.name, ".ExpressionArray.customCDF.txt", sep=""), sep="\t", col.names=F, row.names=F, quote=FALSE)


