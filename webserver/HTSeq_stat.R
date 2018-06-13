

  argv <- commandArgs(TRUE)
  options(stringsAsFactors=F)

# Directory to the htseq result
  htseq.path <- as.character(argv[1])

# Directory to the file which contains list of immune receptor genes
  immuneGene.path <- as.character(argv[2])

# Directory to the final result
  result.path <- as.character(argv[3])
  

  #************************************************************************************
  # Merge multiplt HTSeq results

  htseq_merge <- function(data.path){
      # Arguments:
      #    data.path:path to the individual HTSeq result files  
    
      setwd(data.path)
      options(stringsAsFactors=F)
      files <- list.files()
      files <- files[grep("txt", files)]

      #files <- files[grep("SRX", files)]
      for (file in files) {
          cat(file, "\n")
          temp <- read.table(file, row.names=1, header=F, sep="")
          if (ncol(temp)==0) temp <- read.table(file, row.names=1, header=F, sep="\t")
          if (file==files[1]) {
              counts <- temp
           #    col <- ncol(temp)
           #    row <- nrow(temp)
          } else {
              counts <- cbind(counts, temp[rownames(counts), ])
           #   col <- c(col, ncol(temp))
           #   row <- c(row, nrow(temp))
          }
      }
      colnames(counts) <- gsub("\\.txt", "", files)
      counts

  }

  #************************************************************************************
  # merge the immune receptor gene subtypes

  receptor_merge <- function(expression, gene.list=receptor.ensemble.merge){
      # Arguments:
      #    expression:input expression matrix  
      #    gene.list: list of immune receptor genes

      gene.total <- as.character(unlist(gene.list))
      nonreceptor.expression <- expression[setdiff(rownames(expression), gene.total), ]
      family <- names(gene.list)

      receptor.expression <- c()
      for (gene in family) {
          gene.temp <- as.character(unlist(gene.list[gene]))
          gene.temp <- intersect(gene.temp, rownames(expression))
          cat("The family number of ", gene, " is ", length(gene.temp), "\n")
          if (length(gene.temp) > 1) {
              expression.temp <- apply(expression[gene.temp, ], 2, sum)
          } else {
              expression.temp <- expression[gene.temp, ]
          }
          receptor.expression <- rbind(receptor.expression, expression.temp)          
      }
      rownames(receptor.expression) <- family
      expression <- rbind(nonreceptor.expression, receptor.expression)
      expression
  }

  #************************************************************************************
  # quartile normalization
  quartile <- function(filename, p){
      # Arguments:
      #    filename:input file name  
      #    p: quantile
      #    outname: output file name

      options(stringsAsFactors=F)
      if (is.data.frame(filename)|is.matrix(filename)) {
          data <- filename
      } else {
          if (file.exists(filename)) {
              if (grep("csv", filename)==1)      { data <- read.csv(filename, row.names=1)} 
              else if (grep("txt", filename)==1) { data <- read.csv(filename, row.names=1)}
              else {break}
          }
      }

      n <- ncol(data)
      result <- matrix(nrow=nrow(data), ncol=0)
      for (i in seq(n)) {
          expres <- as.numeric(data[, i])
          expres.min <- min(expres)
          value <- expres[expres != expres.min]
          value <- sort(value, decreasing=T)
          value.quantile <- quantile(value, p)
          scale <- as.numeric(value.quantile)/1000
          expres <- ceiling(expres/scale)
          result <- cbind(result, expres)
      }

      colnames(result) <- colnames(data)
      rownames(result) <- rownames(data)
      result
  }

  #################################################################################

  # 
  counts.raw <- htseq_merge(htseq.path)
  name <- "MouseTissue"
  write.csv(counts.raw, paste(result.path, name, ".Raw.HTSeqData.csv", sep=""), row.names=T)
  n <- nrow(counts.raw)
  index <- grep("ENS", rownames(index))
  counts.raw <- counts.raw[index, ]

  #
  load(immuneGene.path)
  counts.merge <- receptor_merge(counts.raw, gene.list=receptor.ensemble.merge)
  write.csv(counts.merge, paste(result.path, name, ".HTSeqData.immuereceptorMerged.csv", sep=""), row.names=T)

  # normalization
  # The normaliazed data named "counts.quatile" is the input file needed for our web server
  p <- 0.75
  counts.quatile <- quartile(counts.merge, p)
  write.csv(counts.quatile, paste(result.path, name, ".HTSeqData.immuereceptorMerged.quartileNorm.csv", sep=""), row.names=T)
  
