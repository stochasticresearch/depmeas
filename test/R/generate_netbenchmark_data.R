## A script which generates a bunch of data that we can then easily offload to matlab
## for testing hte CIM's ability to infer networks.  

rm(list = ls())
cat("\014")

library(netbenchmark)
library(minet)

netbenchmark_data_gen <- function(datasources.names="all",
                         experiments=150,datasets.num=5,local.noise=20,
                         global.noise=0,noiseType="normal",sym=TRUE,
                         seed=NULL,verbose=TRUE,outputFolder="/home/kiran/ownCloud/PhD/sim_results/netbenchmark/data/")
{
    options(warn=1)
    Fast <- get("Fast", ntb_globals)
    All <- get("All",ntb_globals)
    # set random number generator seed if seed is given
    if (!is.null(seed)) {
        set.seed(seed)
    }else{
        seed <- as.double(Sys.time())
        set.seed(seed)
    }
    
    if(length(datasources.names)==1){
        if (tolower(datasources.names)=="all"){
            datasources.names <- c("rogers1000","syntren1000","syntren300",
                                   "gnw1565","gnw2000")
        }else{
            if(tolower(datasources.names)=="toy"){
                datasources.names <- "toy"
            }
        }
    }
    Availabledata <- get("Availabledata",
                         envir=as.environment("package:grndata"))
    #Availabledata <- eval(parse(text="Availabledata"))
    if(!all(datasources.names %in% Availabledata)){
        stop("The specified datasources are not available")
    }
    seeds <- as.list(round(runif(length(Availabledata),max=1e9)))
    names(seeds) <- Availabledata
    ndata <- length(datasources.names)
    if(!is.na(experiments)){
        for(n in seq_len(ndata)){
            #loading the whole datasource
            datasource <- grndata::getData(datasources.names[n],getNet=FALSE)
            if(experiments*datasets.num>dim(datasource)[1]){
                warning(paste("The specified number of experiments and 
datasets is bigger than the orginal number of experiments in the datasource: 
",datasources.names[n],", sampling with replacement will be used",sep=""))
            } 
        }
    }
    nnets <- datasets.num*ndata
    plots <- list(ndata)
    
    for(n in seq_len(ndata)){ #for each of the datasources
        if(verbose){
            message(paste("Generating DS[",datasources.names[n],"]"))
        }
        aux <- grndata::getData(datasources.names[n])
        datasource <- aux[[1]]
        true.net <- aux[[2]]
        l.seed <- eval(parse(text=paste("seeds$",datasources.names[n])))
        set.seed(l.seed)
        data.list <- datasource.subsample(datasource,
                                          datasets.num=datasets.num,
                                          experiments=experiments,
                                          local.noise=local.noise,
                                          global.noise=global.noise,
                                          noiseType=noiseType)
        
        outputFname = sprintf("%s.Rdata", datasources.names[n])
        fullPathOut = file.path(outputFolder,outputFname)
        #save(list=c("datasource","true.net","data.list"), file = fullPathOut)
        save(list=ls(all.names = TRUE), file = fullPathOut)
    }
}

availableDataSources = c("rogers1000","syntren300","syntren1000","gnw1565","gnw2000")
for(ds in availableDataSources) {
  netbenchmark_data_gen(datasources.names=ds,
                        datasets.num=150,
                        outputFolder="/home/kiran/data/netbenchmark/inputs")
}