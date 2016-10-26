
# Ref
# Gong, T. et al. Optimal deconvolution of transcriptional profiling data
# using quadratic programming with application to complex clinical blood
# samples. PLoS ONE 6, e27156 (2011)


QP <- function(signature, mixture, use.scale='T', QN="T"){
# Solves an lsei inverse problem (Least Squares with Equality and Inequality Constraints)
  library(quadprog)
  library(limSolve)
  library(pcaMethods)
  library(preprocessCore)

# Expression profile format standarlization
  if(rownames(signature)[1] %in% rownames(mixture)){
    data <- data.matrix(mixture)
    signature <- data.matrix(signature)
  }else{
    stop("Input data were in different rowname formats")
  }

  if(max(data) < 50) {data <- 2^data}   
  if(max(signature) < 50) {signature <- 2^signature}  

  if(QN == TRUE){
    tmpc <- colnames(data)
    tmpr <- rownames(data)
    data <- normalize.quantiles(data)
    colnames(data) <- tmpc
    rownames(data) <- tmpr
  } 

  common.name <- intersect(rownames(data), rownames(signature))
  subdata <- data[common.name, ]
  signature <- signature[common.name, ]
 
# number of cell types
  cell.number <- ncol(signature)

## quadratic programming preparation
  if (use.scale) {
     AA <- scale(signature)
  } else {
     AA <- signature
  }
  EE <- rep(1, cell.number)
  FF <- 1
  GG <- diag(nrow=cell.number)
  HH <- rep(0, cell.number)

  out.all <- c()
  sample.number <- ncol(subdata)
  value <- c()
  for (i in seq(sample.number)) {
    BB <- subdata[, i]
    if (use.scale) {
       BB <- scale(BB)
    }  
    out <- lsei(AA, BB, EE, FF, GG, HH)
    u <- sweep(AA, MARGIN=2, out$X, '*')
    v <- apply(u, 1, sum)
    value <- cbind(value, c(v, subdata[, i]))
    rsem <- sqrt(mean((v-BB)^2))
    cat ("The RSEM value is", rsem, "\n")
    corrv <- cor(v, BB)
    cat ("The correlation value is", corrv, "\n")
    out.all <- rbind(out.all, c(out$X * 100, rsem, corrv))
  }
  colnames(out.all) <- c(colnames(AA), "RSEM", "Correlation")
  rownames(out.all) <- colnames(subdata)
  #write.csv(out.all, "QP_Estimated result.csv") 
  out.all
} 
