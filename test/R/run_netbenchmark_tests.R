## R-Script which loads the generated data, and runs it against
## the estimators we care about
rm(list = ls())
cat("\014")

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

MImm <- function(data){ 
  mim <- build.mim(data, estimator = "mi.mm", disc = "equalfreq")
  net <- mrnet(mim)
} 

MIshrink <- function(data){ 
  mim <- build.mim(data, estimator = "mi.shrink", disc = "equalfreq")
  net <- mrnet(mim)
} 

MIsg <- function(data){ 
  mim <- build.mim(data, estimator = "mi.sg", disc = "equalfreq")
  net <- mrnet(mim)
} 

Pearson <- function(data){ 
  mim <- build.mim(data, estimator = "pearson")
  net <- mrnet(mim)
} 

localNoiseVec = c(0,10,20,30,40,50)
globalNoiseVec = c(0,10,20,30,40,50)

for(ds in c("syntren300")) {
    for(gn in globalNoiseVec) {
        for(ln in localNoiseVec) {
            inputFname = sprintf("%s.Rdata", ds)
            subFolder = sprintf('gn_%d_ln_%d',gn,ln)
            fullPathIn = file.path("/data/netbenchmark/inputs",subFolder,inputFname)
            load(file=fullPathIn)
            
            for(i in seq_along(data.list)){
              message('Processing ', ds, '[', i, '/',  length(data.list), ']')
              inputData <- as.data.frame(data.list[[i]])
              #top20.aupr <- netbenchmark.data(methods=c("MIEmpirical", "MImm", "MIshrink", "MIsg", "Kendall", "Pearson"),
              #                                data = inputData,
              #                                true.net=true.net,plot=FALSE,
              #                                verbose=FALSE)
              top20.aupr <- netbenchmark.data(methods=c("Kendall"),
                                              data = inputData,
                                              true.net=true.net,plot=FALSE,
                                              verbose=FALSE)
              # save the output
              outputFname = sprintf("%s_%d.Rdata", ds, i)
              dir.create(file.path("/home/kiran/data/netbenchmark/r_outputs",subFolder), showWarnings = FALSE)
              fullPathOut = file.path("/home/kiran/data/netbenchmark/r_outputs",subFolder,outputFname)
              save(list=c("top20.aupr"), file = fullPathOut)
            }
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
