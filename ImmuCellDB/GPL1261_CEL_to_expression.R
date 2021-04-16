

  #!/usr/bin/Rscript --slave
  cat("********************************************************************************************************************************************\n")
  argv <- commandArgs(TRUE)
  options(stringsAsFactors=F)

  library(mouse4302mmentrezgcdf)
  # mouse4302 mm entrezg cdf version 19.0

  library(affy)
  library(frma)
  library(annotate)
  library(mouse4302mmentrezg.db)

  Data.dir <- as.character(argv[1])
  Result.name <- as.character(argv[2])
  Result.dir <- as.character(argv[3])
  setwd(Data.dir)
  cat(Result.name, "\n")
 
  affydata <- ReadAffy(cdfname = "mouse4302mmentrezgcdf")
  save(affydata, file=paste(Result.dir, Result.name, "_AffyData.RData", sep=""))
  cat("ReadAffy has accompanished!\n")

  # Frma normalization!
  eset <- frma(affydata, normalize="quantile", output.param = "mouse4302mmentrezgfrmavecs")
  cat("FRMA normalization has accompanished!\n")
  save(eset, file = paste(Result.dir, Result.name, "_eset.RData", sep=""))

  ###########################################################
  # Calculation the global normalized unscaled standard error(GNUSE)
  gnuse <- GNUSE(eset, type = "stat")
  colnames(gnuse)<- sub("_.*", "", colnames(gnuse))
  colnames(gnuse) <- sub("\\..*", "", colnames(gnuse))

  #gnuseoutput <- rbind(gnuse, class = rep(celltype, dim(gnuse)[2])) 
  gnusecutoff <- gnuse[1, ] <= 1.25
  gnuse.selected <- gnuse[, gnusecutoff]
  write.csv(t(gnuse.selected), file=paste(Result.dir, Result.name, "_gnuse_quatified.csv", sep=""))  
  save(gnuse, gnuse.selected, file = paste(Result.dir, Result.name, "_GNUSE.RData", sep=""))

  rm(affydata)

  # Obtain gene expression value
  expression <- exprs(eset)
  colnames(expression) <- sub("_.*", "", colnames(expression))
  colnames(expression) <- sub("\\..*", "", colnames(expression))
  expression <- expression[, colnames(gnuse.selected)]
  saveRDS(expression, file=paste(Result.dir, Result.name, "_expression.RDS", sep=""), compress=F)
  cat("Expression data transforming has accompanished!", "\n")

  # Mapping the Probeset into the correspondent EntrezGene
  ID <- featureNames(eset)
  GS <- as.matrix(getSYMBOL(ID, 'mouse4302mmentrezg.db'))    # Mapping the Probeset into the correspondent Gene Symbol
  EG <- as.matrix(getEG(ID, 'mouse4302mmentrezg.db'))   

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

