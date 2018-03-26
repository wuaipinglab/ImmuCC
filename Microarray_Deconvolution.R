
##########################################################################################################################################
###                                                          Description of ImmuCC.R                                            ################
##########################################################################################################################################

  # ImmuCC R script
#       ImmuCC is a tissue deconvolution tool for mouse model that derived from the CIBERSORT method. It was composed of a traing signature matrix specific for mouse and the SVR method called from CIBERSORT.

# 1,    Researchers who are intersted in the application of ImmuCC to estimate the immune composition of mouse tissues, 
#       please cite "Chen,Z.et al.Inference of immune cell composition on the expression profiles of mouse tissue.Sci.Rep.7,40508;doi:10.1038/srep40508(2017)."

# 2,    The core algorithm of ImmuCC is based on a SVR method in CIBERSORT that was developed by Newman et al.
#       Researchers who are intersted in the methodology framework or algorithm, please cite "Newman, A.M.et al. Robust enumeration of cell 
#       subsets from tissue expression profiles. Nature methods 12,453???457,doi:10.1038/nmeth.3337(2015)".

# 3,    To access the CIBERSORT software, please request from https://cibersort.stanford.edu/ and follow their license: http://cibersort.stanford.edu/CIBERSORT_License.txt
#
# 4,    Author: Ziyi Chen, chziy429@163.com, 2017-02-07

##########################################################################################################################################
###                                                          Scripts of ImmuCC                                            ################
##########################################################################################################################################

# Main function of ImmuCC

ImmuCC <- function(path, training_data='srep40508-s1.csv'){

# Function description: 
#     All cel files listed under the specified path can be calculated with function 'ImmuCC'.

# Args:
 #    path: a character denoting the path ReadAffy should look for cel files

 #    training_data: signature matrix for deconvlolution.(The training data srep40508-s1.csv can be 
 #                   downloaded from the supplentary material of Sci.Rep.7,40508;doi:10.1038/srep40508(2017))

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

# Load the function of CIBERSORT
  source("CIBERSORT.R")
# Note: the scirpts of CIBERSORT.R is a method developed by Newman et al.and can be accesssed upon an request from https://cibersort.stanford.edu/

  perm <- 100
  results <- CIBERSORT(training_data, 'mixture.txt', perm)
  results
}
