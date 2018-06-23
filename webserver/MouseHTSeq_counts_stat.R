

  argv <- commandArgs(TRUE)
  options(stringsAsFactors=F)

# Directory to the htseq result
  data.path <- as.character(argv[1])

# path to "receptor.ensemble.merge.RData"
  receptor.ensemble.merge <- as.character(argv[2])

# # Directory to the final normaliazed result
  result.path <- as.character(argv[3])

  #************************************************************************************
  #  Merge the individual htseq result files into one
  htseq_merge <- function(){
      # Arguments:

      options(stringsAsFactors=F)
      files <- list.files()
      files <- files[grep("txt", files)]

      for (file in files) {
          cat(file, "\n")
          temp <- read.table(file, row.names=1, header=F, sep="")
          if (ncol(temp)==0) temp <- read.table(file, row.names=1, header=F, sep="\t")
          if (file==files[1]) {
              counts <- temp
              col <- ncol(temp)
              row <- nrow(temp)
          } else {
              counts <- cbind(counts, temp[rownames(counts), ])
              col <- c(col, ncol(temp))
              row <- c(row, nrow(temp))
          }
      }
      colnames(counts) <- gsub("\\.txt", "", files)
      counts
  }

  #************************************************************************************
  # Merge immune receptor genes into one
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
  # quantile normalization
  quartile <- function(filename, p) {
      # Arguments:
      #    filename:input file name  
      #    p: quantile

      options(stringsAsFactors=F)
      if (is.data.frame(filename)|is.matrix(filename)) {
          data <- filename
      } else {
          if (file.exists(filename)) {
              if (grep("csv", filename)==1)      { data <- read.csv(filename, row.names=1)} 
              else if (grep("txt", filename)==1) { data <- read.table(filename, row.names=1,header=T)}
              else { break }
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

  setwd(data.path)
  if ("receptor.ensemble.merge.RData" %in% list.files()) {
      load(receptor.ensemble.merge)
  } else {
      cat("receptor.ensemble.merge.RData was not existed in the directory, please input the right path or download the file again!", "\n")
  }
  setwd(data.path)

  counts.raw <- htseq_merge()
  # write.csv(counts.raw, paste(name, ".HTSeqData.csv", sep=""), row.names=T, quote=F)

  n <- nrow(counts.raw)
  gene <- grep("ENS", rownames(counts.raw))
  counts.raw <- counts.raw[gene, ]

  counts.merge <- receptor_merge(counts.raw, gene.list=receptor.ensemble.merge)
  # write.csv(counts.merge, paste(name, ".HTSeqData.immuereceptor.csv", sep=""), row.names=T, quote=F)

  p <- 0.75
  counts.quatile <- quartile(counts.merge, p)
  write.csv(counts.quatile, paste(result.path, "/HTSeqData.immuereceptor.quartileNorm.csv", sep=""), row.names=T, quote=F)

