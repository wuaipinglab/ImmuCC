
  argv <- commandArgs(TRUE)
  options(stringsAsFactors=F)

  ## out file name of strand information
  file <- as.character(argv[1])
  outfile <- as.character(argv[2])

  data <- read.delim2(file, header=F)

  if (length(grep("Unknown Data type", data[, 1]))>1){
      cat("There are unknown type files, please check the StrandInfo file!", "\n")
      break
  } else {
      n <- nrow(data)
      samplename <- data[seq(1, n, 5), 1]
      strand1 <- data[grep('\\+\\+,\\-\\-', data[, 1]), 1]
      strand2 <- data[grep("\\+\\-,\\-\\+", data[, 1]), 1]

      strand1 <- unlist(lapply(strand1, function(x){y <- strsplit(x, ": ")[[1]][2]; return (y)}))
      strand2 <- unlist(lapply(strand2, function(x){y <- strsplit(x, ": ")[[1]][2]; return (y)}))

      strand1 <- as.numeric(strand1)
      strand2 <- as.numeric(strand2)

      strandness <- rep("no", length(strand1))
      strandness[strand1>0.8] <- "reverse"
      strandness[strand1>0.8] <- "yes"

      result <- data.frame(samplename=samplename,
                           strand1=strand1,
                           strand2=strand2,
                           strandness=strandness
      )
      write.table(result, outfile, quote=F, row.names=F, col.names=F)

  }
  
