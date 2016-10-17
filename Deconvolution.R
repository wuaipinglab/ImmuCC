

#load essential packages
library(e1071)
library(parallel)
library(preprocessCore)

#Core algorithm
CoreAlg <- function(X, y){
    
    #try different values of nu
    svn_itor <- 3
    
    res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
    }

    if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
    
    nusvm <- rep(0,svn_itor)
    corrv <- rep(0,svn_itor)
    
    #do cibersort
    t <- 1
    while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- cor(k, y)
        t <- t + 1
    }
    
    #pick best model
    rmses <- nusvm
    mn <- which.min(rmses)
    model <- out[[mn]]
    
    #get and normalize coefficients
    q <- t(model$coefs) %*% model$SV
    q[which(q<0)]<-0
    w <- (q/sum(q))
    
    mix_rmse <- rmses[mn]
    mix_r <- corrv[mn]

    newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#do permutations
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#main function
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE, samplename){
  #read in data
  X <- read.csv(sig_matrix, header=T, row.names=1, check.names=F)
  Y <- read.csv(mixture_file, header=T, row.names=1,check.names=F)

  Y<-Y[rownames(X)[rownames(X) %in% rownames(Y)],]
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  } 
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
  #print(nulldist)
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor], w*100, pval, mix_r, mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  write.table(rbind(header, output), file=paste(samplename, "_.txt", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
    obj <- output[, -1]
   # obj <- data.matrix(obj)
   # rownames(obj) <- colnames(Y)
   # colnames(obj) <- header
   # write.csv(obj, file=paste(samplename,"_CIBERSORT-Results.csv"))
    return(rbind(header, output))
}


# Transform the raw microarray CEL files into expressionSet
cel_expression <- function(){
    library(affy)
    library(frma)
    library(org.Mm.eg.db)
    library(mouse4302frmavecs)
    library(mouse4302mmentrezgcdf)
    library(mouse4302mmentrezg.db)
    # read all CEL files in current directory
    affydata <- ReadAffy(cdfname = "mouse4302mmentrezgcdf")
    # FRMA(Background correction normalization and summarization)
    eset <- frma(affydata)
    # Construct mapping matrix
    ID <- featureNames(eset)
    GS <- as.matrix(getSYMBOL(ID, 'mouse4302mmentrezg.db'))    # Mapping the Probeset into the correspondent Gene Symbol
    cols =  c("ProbeSet", "GeneSymbol")
    mapping <- cbind(ID, GS)
    mapping <- mapping[which(mapping[, 2] != "NA"), ]            # remove NAs
    colnames(mapping) <- cols
    rownames(mapping) <- mapping[, 1]

    # Calculation the global normalized unscaled standard error(GNUSE)
    gnuse <- GNUSE(eset, type = "stat")
    gnuseoutput <- rbind(gnuse, class = rep(celltype, dim(gnuse)[2])) 
    gnusecutoff <- gnuse[1, ] <= 1.25  

    # Obtain gene expression profile
    express <- exprs(eset)
    express <- express[rownames(mapping), gnusecutoff]
    colnames(express) <- sub("_.*", "", colnames(express))
    colnames(express) <- sub("\\..*", "", colnames(express))
    rownames(express) <- mapping[rownames(express), 2]
    save(express, file="Expression.RData")

    # save the result
    cat("Expression data transforming has accompanished!", "\n")
    write.csv(express, file="ExpressionData.csv", col.names=T, row.names=F)
}

dir <- ""
# set work directory
setwd(dir)
# Convert raw CEL microarray files into expressionSet
cel_expression()
output_name <- ""
options(stringsAsFactors=F)
## Main function
main <- function(mixture_name="ExpressionData.csv", sig_name="signature_matrix(GS).csv", output_name){ 
# Argument:
#     mixture_name: input mixture sample name
#     sig_name: input signature matrix name
#     output_name: output file name

  # load in the expression data!
  expression <- read.csv(mixture_name, header=T, row.names=1) 
  # load in the signature matrix!
  result <- CIBERSORT(sig_matrix=sig, mixture_file = expression, perm=100, QN=T, output_name)
  
}


















