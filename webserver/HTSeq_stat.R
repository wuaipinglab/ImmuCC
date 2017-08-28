

  argv <- commandArgs(TRUE)
  options(stringsAsFactors=F)

  #************************************************************************************
  htseq_merge <- function(){
      # Arguments:
      #    path:path to the input files  

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
           #   col <- ncol(temp)
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
          cat("The family number of ", gene, " is ",length(gene.temp), "\n")
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
          value <- expres[expres!=0]
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
  setwd(as.character(argv[1]))

  counts.raw <- htseq_merge()
  name <- "MouseTissue"
  write.csv(counts.raw, paste(name, ".HTSeqData.csv", sep=""), row.names=T, quote=F)

  n <- nrow(counts.raw)
  index <- grep("ENS", rownames(index))
  counts.raw <- counts.raw[index, ]
  counts.raw <- counts.raw[1:(n-5), ]

  load("./receptor.ensemble.merge.RData")
  counts.merge <- receptor_merge(counts.raw, gene.list=receptor.ensemble.merge)
  write.csv(counts.merge, paste(name, ".HTSeqData.immuereceptor.csv", sep=""), row.names=T, quote=F)

  p <- 0.75
  counts.quatile <- quartile(counts.merge, p)
  write.csv(counts.quatile, paste(name, ".HTSeqData.immuereceptor.quartileNorm.csv", sep=""), row.names=T, quote=F)
  
