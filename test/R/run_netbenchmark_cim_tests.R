## R-Script which loads the generated data, and runs it against
## the estimators we care about
rm(list = ls())
cat("\014")

library(netbenchmark)
library(minet)
readMat <- R.matlab::readMat

#' Benchmarking of several network inference algorithms for your own data
#' modified to read in a file where the MIM was computed on a different machine/processor
netbenchmark_custom.data <- function(methods="all.fast",data=NULL,true.net=NULL,
                              eval="AUPR",no.topedges=20,sym=TRUE,plot=FALSE,
                              verbose=TRUE)
{
  options(warn=1)
  Fast <- get("Fast", ntb_globals)
  All <- get("All",ntb_globals)
  if(all("all.fast" %in% tolower(methods))) {
    methods <- c(Fast,methods[tolower(methods)!="all.fast"])
  }else if(all("all" %in% tolower(methods))) {
    methods <- c(All,methods[tolower(methods)!="all"])
  }
  if(is.null(data)){
    stop("You should provide the data.")
  }
  if(is.null(true.net)){
    stop("You should provide the underlying network.")
  }
  # if(!is.data.frame(data)){
  #     stop("The provided data shoud be a data.frame.")
  # }
  # if(!is.matrix(true.net)){
  #     stop("The provided true.net shoud be a adjacency matrix.")
  # }
  # if(dim(true.net)[1]!=dim(true.net)[2]){
  #     stop("The provided true.net shoud be a adjacency matrix.")
  # }
  # if(all(dim(data)[2]!=dim(true.net)[1])){
  #     stop("The provided data shoud contain variables in columuns")
  # }
  # if(all(colnames(data)!=colnames(true.net)[1])){
  #     stop("The provided data shoud contain variables in columuns")
  # }
  nnets <- 1
  ndata <- 1
  nmeths <- length(methods)
  results.table <-  as.data.frame(matrix(0,nrow=nnets,ncol=nmeths+1))
  nlinks.table <- matrix(0,nrow=nnets,ncol=nmeths+1)
  time.table <-  as.data.frame(matrix(0,nrow=nnets,ncol=nmeths))
  names <- as.character(methods)
  colnames(results.table) <- c(names,"rand")
  rownames(results.table) <- as.character(1:nnets)
  colnames(time.table) <- c(names)
  rownames(time.table) <- as.character(1:nnets)
  plots <- list(ndata)
  #########################
  ngenes <- dim(true.net)[1]
  npos <- sum(true.net) #number of true links in the network
  nlinks <- ngenes^2-ngenes #number of posible links in the network
  if(sym){
    nlinks <- nlinks/2
  }
  no.edges <- round(nlinks*no.topedges/100)
  nneg <- nlinks-npos #number of false links in the network
  conf.mat <-  matrix(0,nrow=nmeths,ncol=4)
  colnames(conf.mat) <- c("TP","FP","TN","FN")
  rownames(conf.mat) <- names
  best <- -1
  if(plot){
    col <- rainbow(nmeths)
  }
  for(j in seq_len(nmeths)){
    if(verbose){
      message(names[j])
    }    
    ptm <- proc.time()
    net <- do.call(names[j],list(data, true.net))
    t <- proc.time() - ptm
    if(sum(is.na(net)>0)){
      net[is.na(net)] <- 0
    }
    
    r <- evaluate(net,true.net,extend=no.edges,sym=sym)
    time.table[1,j]=t[[3]]
    conf.mat[j,] <- r[no.edges,]
    if( tolower(eval)=="no.truepos"){
      results.table[1,j]=mean(r[1:no.edges,"TP"])
    }else if (tolower(eval)== "aupr"){
      results.table[1,j]=aupr(r,no.edges)
    }else if (tolower(eval)== "auroc"){
      results.table[1,j]=auroc(r,no.edges)
    }else stop("unknown evaluation metric")
    if(results.table[1,j]>best){
      best <- results.table[1,j]
      best.net <- net
      M <- j
    }
    if(plot){
      col<-rainbow(nmeths)
      if(j==1){
        pr.plot(r[1:no.edges,],col=col[j],lwd=2)
      }else{
        pr.plot(r[1:no.edges,],col=col[j],lwd=2,device=2)
      }
    }
    table <- as.data.frame(r[1:no.edges,])
    names(table) <- sapply(names(table),tolower)
    plots[[j]] <-  minet::pr(table)
  }
  if(M!=which.max(results.table[1,1:(nmeths)])){
    stop("error")
  }
  rand.net <- matrix(runif(ngenes^2),ngenes,ngenes)
  diag(rand.net) <- 0
  colnames(rand.net) <- colnames(true.net)
  rownames(rand.net) <- colnames(true.net)
  r <- evaluate(rand.net,true.net,extend=no.edges,sym=sym)
  if( tolower(eval)== "no.truepos"){
    results.table[1,nmeths+1]=mean(r[1:no.edges,"TP"])
  }else if (tolower(eval)== "aupr"){
    results.table[1,nmeths+1]=aupr(r,no.edges)
  }else if (tolower(eval)== "auroc"){
    results.table[1,nmeths+1]=auroc(r,no.edges)
  }
  L <- list(results.table,time.table,plots)
  names(L) <- c(paste(eval,"top",as.character(no.topedges),"%",sep=""),
                "CpuTime","PRcurves")
  return(L)
}

CIM <- function(ds, true.net){ 
  
  # load the MIM matrix that was processed in Matlab
  mim <- readMat(ds[1])
  # perform pre-processing on the R matrix
  mim <- (mim$R)^2
  diag(mim) <- 0
  maxi<-0.999999
  mim[which(mim>maxi)]<-maxi
  mim <- -0.5*log(1-mim)
  mim[mim<0]<-0
  
  # assign the names so that evaluate works properly
  colnames(mim) <- colnames(true.net)
  rownames(mim) <- rownames(true.net)
  
  # construct the mrnet network
  net <- mrnet(mim)
} 

MatlabMI <- function(ds, true.net) {
  # load the MIM matrix that was processed in Matlab
  mim <- readMat(ds[1])$R
  
  # assign the names so that evaluate works properly
  colnames(mim) <- colnames(true.net)
  rownames(mim) <- rownames(true.net)
  
  # construct the mrnet network
  net <- mrnet(mim)
}

availableDataSources = c("syntren300","rogers1000","syntren1000","gnw1565","gnw2000")

if(Sys.info()["sysname"]=="Darwin") {  # mac
  inputRepo    = "/Users/kiran/data/netbenchmark/inputs"
  rResultsRepo = "/Users/kiran/data/netbenchmark/r_outputs"
  matlabResultsRepo = "/Users/kiran/data/netbenchmark/matlab_outputs"
} else {  # unix
  inputRepo    = "/home/kiran/data/netbenchmark/inputs"
  rResultsRepo = "/home/kiran/data/netbenchmark/r_outputs"
  matlabResultsRepo = "/home/kiran/data/netbenchmark/matlab_outputs"
}

for(ds in availableDataSources) {
  # load the datasource to get the true network
  inputFname = file.path(inputRepo,sprintf("%s.Rdata",ds))
  load(inputFname)
  for(i in seq_along(data.list)){
    cimInputFname   = sprintf("%s_%d_cim_output.mat", ds, i)
    knn1InputFname  = sprintf("%s_%d_knn1_output.mat", ds, i)
    knn6InputFname  = sprintf("%s_%d_knn6_output.mat", ds, i)
    knn20InputFname = sprintf("%s_%d_knn20_output.mat", ds, i)
    vmeInputFname   = sprintf("%s_%d_vme_output.mat", ds, i)
    apInputFname    = sprintf("%s_%d_ap_output.mat", ds, i)
    # TODO: do TAU if you see any anomalous results ...

    cim_top20.aupr   <- netbenchmark_custom.data(methods=c("CIM"),data = file.path(matlabResultsRepo,cimInputFname),
                                    true.net=true.net,plot=FALSE,verbose=FALSE)
    knn1_top20.aupr  <- netbenchmark_custom.data(methods=c("MatlabMI"),data = file.path(matlabResultsRepo,knn1InputFname),
                                    true.net=true.net,plot=FALSE,verbose=FALSE)
    knn6_top20.aupr  <- netbenchmark_custom.data(methods=c("MatlabMI"),data = file.path(matlabResultsRepo,knn6InputFname),
                                    true.net=true.net,plot=FALSE,verbose=FALSE)
    knn20_top20.aupr <- netbenchmark_custom.data(methods=c("MatlabMI"),data = file.path(matlabResultsRepo,knn20InputFname),
                                    true.net=true.net,plot=FALSE,verbose=FALSE)
    vme_top20.aupr   <- netbenchmark_custom.data(methods=c("MatlabMI"),data = file.path(matlabResultsRepo,vmeInputFname),
                                    true.net=true.net,plot=FALSE,verbose=FALSE)
    ap_top20.aupr    <- netbenchmark_custom.data(methods=c("MatlabMI"),data = file.path(matlabResultsRepo,apInputFname),
                                    true.net=true.net,plot=FALSE,verbose=FALSE)
    
    print(cim_top20.aupr$`AUPRtop20%`)
    #print(knn1_top20.aupr$`AUPRtop20%`)

    # save the output
    outputFname = sprintf("%s_%d_matlab.Rdata", ds, i)
    fullPathOut = file.path(rResultsRepo,outputFname)
    save(list=c("cim_top20.aupr","knn1_top20.aupr","knn6_top20.aupr","knn20_top20.aupr",
                "vme_top20.aupr","ap_top20.aupr"),file = fullPathOut)
    #save(list=c("cim_top20.aupr"),file = fullPathOut)
  }
}

# syntrenIdx = 1
# 
# fullPathIn = file.path(matlabResultsRepo,sprintf("syntren300_%d_output.mat",syntrenIdx))
# z <- readMat(fullPathIn)
# mimCim <- (z$R)^2
# diag(mimCim) <- 0
# maxi<-0.999999
# mimCim[which(mimCim>maxi)]<-maxi
# mimCim <--0.5*log(1-mimCim)
# mimCim[mimCim<0]<-0
# # construct the mrnet network
# netCim <- mrnet(mimCim)
# 
# # load the "true" network
# inputFname = file.path(inputRepo,"syntren300.Rdata")
# load(inputFname)
# 
# # assign col/row names
# colnames(netCim) <- colnames(true.net)
# rownames(netCim) <- rownames(true.net)
# 
# # make sure no N/A in the networks
# if(sum(is.na(netCim)>0)){
#   netCim[is.na(netCim)] <- 0
# }
# 
# ngenes <- dim(true.net)[1]
# npos <- sum(true.net) #number of true links in the network
# nlinks <- ngenes^2-ngenes #number of posible links in the network
# no.topedges <- 20
# if(sym){
#   nlinks <- nlinks/2
# }
# no.edges <- round(nlinks*no.topedges/100)
# 
# rCim  <- evaluate(netCim,true.net,extend=no.edges,sym=TRUE)
# # cimCalcLen = min(length(rTau[,"TP"]),length(rCim[,"TP"]))
# # print(sprintf("(tau-cim) = %d %d %d %d",
# #       max(rTau[1:cimCalcLen,"TP"]-rCim[1:cimCalcLen,"TP"]),
# #       max(rTau[1:cimCalcLen,"FP"]-rCim[1:cimCalcLen,"FP"]),
# #       max(rTau[1:cimCalcLen,"TN"]-rCim[1:cimCalcLen,"TN"]),
# #       max(rTau[1:cimCalcLen,"FN"]-rCim[1:cimCalcLen,"FN"])))
# 
# # fullPathIn contains the path to the Matlab Syntren outputs
# cim_top20.aupr <- netbenchmark_custom.data(methods=c("CIM"),data = fullPathIn,
#                                            true.net=true.net,plot=FALSE,verbose=FALSE)
# print(cim_top20.aupr$`AUPRtop20%`)