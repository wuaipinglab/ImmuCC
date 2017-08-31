 
  options(stringsAsFactors=F)
########################################################################################################
##                                    Affy pre-procession                           ################
########################################################################################################

  library(affy)
  library(frma)
  library(mouse4302mmentrezgcdf)
  library(mouse4302frmavecs)
  library(preprocessCore)

# Read all cel files under path with a custom cdf "mouse4302mmentrezcdf"
  affydata <- ReadAffy(celfile.path=path, cdf="mouse4302mmentrezcdf")

# Preprocessing with frma
  eset <- frma(affydata)

# Output the expression value of samples profiled on array
  ematrix <- exprs(eset)
  write.table(ematrix, "mixture.txt",row.names=F, col.names=F)
  write.csv(ematrix, "mixture.csv")


######################################################################################################################
###                              Agilent pre-procession                           ################
######################################################################################################################
platform <- ""
outname <- ""
#Read the raw txt files
library(limma)

# filename: vector of filenames to be read
filename <- list.files()[grep("*.txt.gz", list.files())]
x <- read.maimages(filename, source="agilent", green.only=TRUE)

# Background correct and normalize
y <- backgroundCorrect(x, method="normexp")
y <- normalizeBetweenArrays(y, method="quantile")
expression <-y$E
rownames(expression) <- y$genes$ProbeName
gene <- y$genes$GeneName
# y$E            Gene expression value
# y$targets      file name
# y$genes        Gene and probe names

# mapping information
mapping <- matrix(c(y$genes$Row, y$genes$Col, y$genes$Start, y$genes$Sequence, y$genes$ProbeUID, 
                  y$genes$ControlType, y$genes$ProbeName, y$genes$GeneName, y$genes$SystematicName, y$genes$Description), ncol=10)
save(mapping, file=paste(platform, "_mapping.RData", sep=""))

# mapping information + mean expression value of probes
expressionStat <- cbind(mapping, rowMeans(y$E))

# Unique gene names
UniGene <- unique(y$genes$GeneName)

# Gene deduplicated
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
write.csv(select, paste(platform, "SelectedMapping.csv", sep=""), row.names=F)

rownames(select) <- select[, 7]
expression <- expression[select[, 7], ]
rownames(expression) <- select[rownames(expression), 8]
#save(expression, file=paste(platform, "FilterExpressionSet.RData", sep=""))

write.table(expression, file=paste(outname, "_ExpressionMatrix.txt", sep=""),row.names=F, col.names=F)
write.csv(expression, file=paste(outname, "_ExpressionMatrix.csv", sep=""))
