 
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
