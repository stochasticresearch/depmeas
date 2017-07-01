## R-Script which loads the generated data, and runs it against
## the estimators we care about
rm(list = ls())
cat("\014")

RUN_CIM = 0;

library(netbenchmark)
library(minet)

Kendall <- function(data){ 
  mim <- build.mim(data, estimator = "kendall")
  net <- mrnet(mim)
} 

MIEmpirical <- function(data){ 
  mim <- build.mim(data, estimator = "mi.empirical", disc = "equalfreq")
  net <- mrnet(mim)
} 

MIsg <- function(data){ 
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
    top20.aupr <- netbenchmark.data(methods=c("Kendall", "MIEmpirical", "MIsg"),data = inputData,
                                    true.net=true.net,plot=FALSE)
    # save the output
    outputFname = sprintf("%s.Rdata", ds)
    fullPathOut = file.path(outputFolder,outputFname)
    save(list=c("top20.aupr"), file = fullPathOut)
  }
}

if(RUN_CIM) {
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
}

###############################################################################
# load the data
#fname = '/home/kiran/ownCloud/PhD/sim_results/netbenchmark/data/syntren300.Rdata'
#load(file=fname)
#inputData <- as.data.frame(data.list[[2]])
#top20.aupr <- netbenchmark.data(methods=c("MIEmpirical", "Kendall"),data = inputData,
#                                true.net=true.net,plot=FALSE)
#top20.aupr <- netbenchmark(datasources.names=c("syntren300"), 
#                     methods=c("MIEmpirical", "Pearson"),
#                     datasets.num=150,
#                     verbose=FALSE,
#                     seed = 123)
#print(top20.aupr[[1]])
###############################################################################