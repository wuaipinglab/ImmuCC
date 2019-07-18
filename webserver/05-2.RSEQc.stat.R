
  argv <- commandArgs(TRUE)
  options(stringsAsFactors=F)

  ## out file name of strand information
  file <- argv[1]
  outfile <- argv[2]

  data <- read.delim2(file, header=F)

  if (length(grep("Unknown Data type", data[, 1]))>1){
      cat("There are unknown type files, please check the StrandInfo file!", "\n")
      break
  } else {
      index1 <- grep('This is', data[, 1])
      index2 <- grep('Fraction', data[, 1])
      samplename <- data[-c(index1, index2), 1]
      sample.index <- which(data[, 1] %in% samplename)
      inform.range <- c(sample.index[-1], nrow(data)+1)-c(sample.index)
      samplename <- samplename[inform.range ==5]
      
      strand1 <- data[grep('1\\+\\+,1\\-\\-', data[, 1]), 1]
      strand2 <- data[grep("1\\+\\-,1\\-\\+", data[, 1]), 1]

      strand1 <- unlist(lapply(strand1, function(x){y <- strsplit(x, ": ")[[1]][2]; return (y)}))
      strand2 <- unlist(lapply(strand2, function(x){y <- strsplit(x, ": ")[[1]][2]; return (y)}))

      strand1 <- as.numeric(strand1)
      strand2 <- as.numeric(strand2)
    
      # Samples of low quality. Strandness information is ambiguous
      note <- rep("no", length(strand1))
      note[strand1<0.7 & strand1>=0.6] <- "low Quality"
    
      # Strandness parameter for HTSeq
      strandness1 <- rep("no", length(strand1))
      strandness1[strand2>0.7] <- "reverse"
      strandness1[strand1>0.7] <- "yes"
    
      # Strandness parameter for RSEM
      strandness2 <- rep("none", length(strand1))
      strandness2[strand2>0.7] <- "reverse"
      strandness2[strand1>0.7] <- "forward"

      result <- data.frame(samplename=samplename,
                           strand1=strand1,
                           strand2=strand2,
                           strandness1=strandness1,
                           strandness2=strandness2,
                           Notation=note
      )
      write.table(result, outfile, quote=F, row.names=F, col.names=F)

  }
  
