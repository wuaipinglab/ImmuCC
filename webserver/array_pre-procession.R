 
  options(stringsAsFactors=F)
########################################################################################################
##                                      Affy pre-procession                             ################
########################################################################################################

  library(affy)
  library(frma)
  library(mouse4302mmentrezgcdf)
  library(mouse4302frmavecs)
  library(preprocessCore)

# Directory to the cel files
  cel.path <- ""

# Directory to the result path
  result.path <- ""

# Read all cel files with a custom cdf "mouse4302mmentrezcdf"
  affydata <- ReadAffy(celfile.path=cel.path, cdf="mouse4302mmentrezcdf")

# Preprocessing with frma
  eset <- frma(affydata)

# Get the expression profile of all samples
  ematrix <- exprs(eset)

# svae the expression profile
  write.table(ematrix, paste(result.path, "mixture.txt", sep=""),row.names=F, col.names=F)
  write.csv(ematrix, paste(result.path, "mixture.csv", sep=""))


######################################################################################################################
########                                         Agilent pre-procession                                       ########
######################################################################################################################

options(stringsAsFactors=F)
library(limma)

# Name of the platform
  platform <- ""

# Prefix of the result file name
  outname <- ""

# Directory to the result path
  result.path <- ""

# filename: vector of filenames to be read
  filename <- list.files()[grep("*.txt.gz", list.files())]

#Read the raw txt files
  x <- read.maimages(filename, source="agilent", green.only=TRUE)

# Background correction and normalization
  y <- backgroundCorrect(x, method="normexp")
  y <- normalizeBetweenArrays(y, method="quantile")
  expression <-y$E
  rownames(expression) <- y$genes$ProbeName
  gene <- y$genes$GeneName
# y$E            Gene expression value
# y$targets      file name
# y$genes        Gene and probe names

# Extract the mapping information
  mapping <- matrix(c(y$genes$Row, y$genes$Col, y$genes$Start, y$genes$Sequence, y$genes$ProbeUID, 
                    y$genes$ControlType, y$genes$ProbeName, y$genes$GeneName, y$genes$SystematicName, y$genes$Description), ncol=10)
  save(mapping, file=paste(result.path, platform, "_mapping.RData", sep=""))

# mapping information + probe mean expression values
  expressionStat <- cbind(mapping, rowMeans(y$E))

# Get the unique gene names
  UniGene <- unique(y$genes$GeneName)

# Remove the duplicated genes
  select <- matrix(nrow=0, ncol=11)
  for (gene in UniGene) {
      cat("Curent Gene is:", gene , "\n")
      temp <- expressionStat[expressionStat[, 8]==gene, ]
      temp <- matrix(temp, ncol=11, byrow=F)
      if (dim(temp)[1]==1) {
          select <- rbind(select, temp)
      } else {
          temp <- temp[order(temp[, 11], decreasing=T), ]
          select <- rbind(select, temp[1, ])
      }    
  }

# result saving
  colnames(select) <- c(names(y$genes), "MeanValue")
  write.csv(select, paste(result.path, platform, "SelectedMapping.csv", sep=""), row.names=F)

  rownames(select) <- select[, 7]
  expression <- expression[select[, 7], ]
  rownames(expression) <- select[rownames(expression), 8]

# svae the expression profile
  colnames(expression) <- gsub("\\..*", "", colnames(expression))
  write.table(expression, file=paste(result.path, outname, "_ExpressionMatrix.txt", sep=""),row.names=F, col.names=F)
  write.csv(expression, file=paste(result.path, outname, "_ExpressionMatrix.csv", sep=""))
