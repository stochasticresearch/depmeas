## R-Script which loads the generated data, and runs it against
## the estimators we care about
rm(list = ls())
cat("\014")

library(netbenchmark)
library(minet)

CIM <- function(data){ 
  mim <- build.mim(data, estimator = "mi.sg", disc = "equalfreq")
  net <- mrnet(mim)
} 

availableDataSources = c("syntren300","rogers1000","syntren1000","gnw1565","gnw2000")
dataRepo = "/home/kiran/data/netbenchmark/inputs/"
outputFolder = "/home/kiran/data/netbenchmark/r_outputs"

for(ds in availableDataSources) {
  inputFname = sprintf("%s.Rdata", ds)
  fullPathIn = file.path(dataRepo,inputFname)
  load(file=fullPathIn)
  
  for(i in seq_along(data.list)){
    inputData <- as.data.frame(data.list[[i]])
    top20.aupr <- netbenchmark.data(methods=c("CIM"),data = inputData,
                                    true.net=true.net,plot=FALSE)
    # save the output
    outputFname = sprintf("%s_CIM.Rdata", ds)
    fullPathOut = file.path(outputFolder,outputFname)
    save(list=c("top20.aupr"), file = fullPathOut)
  }
}