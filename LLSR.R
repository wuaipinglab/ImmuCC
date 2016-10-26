

# Ref
# Abbas, A.R., Wolslegel, K., Seshasayee, D., Modrusan, Z. & Clark, H.F.
# Deconvolution of blood microarray data identifies cellular activation
# patterns in systemic lupus erythematosus. PLoS one 4, e6098
# (2009)

LLSR <- function(XX, YY, w=NA, QN="T"){
# Argues:
# XX: signature matrix
# YY: sample expression profile

  library(preprocessCore)
  # Expression profile format standarlization
  XX <- data.matrix(XX)
  YY <- data.matrix(YY)

  if(max(XX) < 50) {XX <- 2^XX}   
  if(max(YY) < 50) {YY <- 2^YY} 
  
  if(QN == T){
    tmpc <- colnames(YY)
    tmpr <- rownames(YY)
    data <- normalize.quantiles(YY)
    colnames(YY) <- tmpc
    rownames(YY) <- tmpr
  } 

  name <- intersect(rownames(XX), rownames(YY))
  XX <- XX[name, ]
  YY <- YY[name, ]

  N1 <- ncol(YY)
  fraction <- c()
  for(i in seq(N1)){
    y=YY[, i]
    if(is.na(w[1])) tmp=lsfit(XX, y, intercept=F) else tmp=lsfit(XX, y, w, intercept=F)
    tmp.beta=tmp$coefficients
    tmp.beta[which(tmp.beta <0)] <- 0
    cat(length(tmp.beta), "\n")
    u <- sweep(XX, MARGIN=2, tmp.beta, '*')
    v <- apply(u, 1, sum)
    corrv <- cor(v, YY[, i])
    fraction <- rbind(fraction, c(tmp.beta * 100, corrv))
  }
  colnames(fraction) <- c(colnames(XX), "Correlation")
  rownames(fraction) <- colnames(YY)
  fraction
}
