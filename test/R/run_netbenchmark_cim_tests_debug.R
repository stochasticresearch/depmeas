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

evaluate_custom <- function(inf.net, true.net,sym=TRUE,extend=0) 
{
  if((diff(dim(inf.net))!=0) && dim(inf.net)[2]==3){
    ngenes <- max(as.double(inf.net[,1:2]))
    if((diff(dim(true.net))!=0) && dim(true.net)[2]==3){
      ngenes <-max(as.double(true.net[,1:2]),ngenes)
      true.net <- Matrix::sparseMatrix(i=true.net[,1],
                                       j=true.net[,2],x=true.net[,3],
                                       dims =c(ngenes,ngenes))
      true.net<-as.matrix(true.net)
    }else{
      ngenes <- dim(true.net)[2]
    }
    inf.net <- Matrix::sparseMatrix(i=inf.net[,1],
                                    j=inf.net[,2],x=inf.net[,3],
                                    dims =c(ngenes,ngenes))
    inf.net<-as.matrix(inf.net)
  }else{
    ngenes <- dim(true.net)[2]
  }
  if(sym){
    if(extend>(ngenes*(ngenes+1)/2)){
      extend <- (ngenes*(ngenes+1)/2)
    }
    true.net <- pmax(true.net,t(true.net))
    inf.net <- pmax(inf.net,t(inf.net))
    if(is.null(colnames(true.net))){
      colnames(true.net) <- as.character(1:ngenes)
      rownames(true.net) <- as.character(1:ngenes)
      colnames(inf.net) <- as.character(1:ngenes)
      rownames(inf.net) <- as.character(1:ngenes)
    }else if(is.null(colnames(inf.net))){
      colnames(inf.net) <- colnames(true.net)
      rownames(inf.net) <- rownames(true.net)
    }
    
    true.net[lower.tri(true.net,diag=TRUE)] <- 0
    inf.net[lower.tri(inf.net,diag=TRUE)] <- 0
    if(sum(inf.net!=0)>0){
      PredEdgeList <- .Adj2Edgelist(inf.net)
    }
    if(sum(inf.net!=0)<extend){
      inv.net <- inf.net*0
      inv.net[idxInvert(inf.net,which(inf.net!=0))] <- 1
      diag(inv.net) <- 0
      inv.net[lower.tri(inv.net)] <- 0
      E <- .Adj2Edgelist(inv.net)
      idx <- dim(E)[1]
      E <- E[sample(idx),]
      if(sum(inf.net!=0)==0){
        print(">> 5a")
        PredEdgeList <- E[1:(extend-sum(inf.net!=0)+1),]
      }else{
        print(">> 5b")
        PredEdgeList <- rbind(PredEdgeList,
                              E[1:(extend-sum(inf.net!=0)),])
      }
    }
    GSEdgeList <- .Adj2Edgelist(true.net)
    r <- rate(PredEdgeList, GSEdgeList, ngenes, 2)
  }else{
    if(extend>(ngenes^2-ngenes)){
      extend <- (ngenes^2-ngenes)
    }
    if(is.null(colnames(true.net))){
      colnames(true.net) <- as.character(1:ngenes)
      rownames(true.net) <- as.character(1:ngenes)
      colnames(inf.net) <- as.character(1:ngenes)
      rownames(inf.net) <- as.character(1:ngenes)
    }else if(is.null(colnames(inf.net))){
      colnames(inf.net) <- colnames(true.net)
      rownames(inf.net) <- rownames(true.net)
    }
    if(sum(inf.net!=0)>0){
      PredEdgeList <- .Adj2Edgelist(inf.net)
    }
    if(sum(inf.net!=0)<extend){
      inv.net <- inf.net*0
      inv.net[idxInvert(inf.net,which(inf.net!=0))] <- 1
      diag(inv.net) <- 0
      E <- .Adj2Edgelist(inv.net)
      idx <- dim(E)[1]
      E <- E[sample(idx),]
      if(sum(inf.net!=0)==0){
        PredEdgeList <- E[1:(extend-sum(inf.net!=0)),]
      }else{
        PredEdgeList <- rbind(PredEdgeList,
                              E[1:(extend-sum(inf.net!=0)),])
      }
    }
    GSEdgeList <- .Adj2Edgelist(true.net)
    r <- rate(PredEdgeList, GSEdgeList, ngenes, 1)
  }
  colnames(r) <- c("TP","FP","TN","FN")
  return(r)
}

.Adj2Edgelist <- function(Adjmat){
  idx <- which(Adjmat!=0, arr.ind=TRUE)
  r=idx[,1]; c=idx[,2];
  E <- as.character(rep(0,dim(idx)[1]*3))
  l <- dim(idx)[1]
  dim(E) <- c(l,3)
  names <- colnames(Adjmat)
  A <- rep(0,l)
  for(i in seq_len(l)){
    E[i,1] <- names[r[i]]
    E[i,2] <- names[c[i]]
    A[i] <- Adjmat[r[i],c[i]]
  }
  a <- sort(A,decreasing = TRUE,index.return = TRUE)
  E <- E[a$ix,]
  dim(E) <- c(l,3)
  E[,3] <- as.character(a$x)
  return(E)
}


Kendall <- function(data){ 
  mim <- build.mim(data, estimator = "kendall")
  net <- mrnet(mim)
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

# for(ds in availableDataSources) {
#   # load the datasource to get the true network
#   inputFname = file.path(inputRepo,sprintf("%s.Rdata",ds))
#   load(inputFname)
#   #for(i in seq_along(data.list)){
#   for(i in seq_along(1)){
#     inputFname = sprintf("%s_%d_output.mat", ds, i)
#     fullPathIn = file.path(matlabResultsRepo,inputFname)
# 
#     #cim_top20.aupr <- netbenchmark_custom.data(methods=c("CIM"),data = fullPathIn,
#     #                                true.net=true.net,plot=FALSE,verbose=FALSE)
#     #print(cim_top20.aupr$`AUPRtop20%`)
#     
#     net <- CIM(c(fullPathIn))
#     image(net)
#     
#     # save the output
#     outputFname = sprintf("%s_%d_CIM.Rdata", ds, i)
#     fullPathOut = file.path(rResultsRepo,outputFname)
#     save(list=c("cim_top20.aupr"), file = fullPathOut)
#   }
# }

syntrenIdx = 1

fullPathIn = file.path(matlabResultsRepo,sprintf("syntren300_%d_output.mat",syntrenIdx))
z <- readMat(fullPathIn)
mimCim <- (z$R)^2
diag(mimCim) <- 0
maxi<-0.999999
mimCim[which(mimCim>maxi)]<-maxi
mimCim <--0.5*log(1-mimCim)
mimCim[mimCim<0]<-0
# construct the mrnet network
netCim <- mrnet(mimCim)

# load the "true" network
inputFname = file.path(inputRepo,"syntren300.Rdata")
load(inputFname)

# try to run against kendalls tau
ds <- data.list[[syntrenIdx]]
mimTau<-cor(ds,method="kendall",use="complete.obs")
mimTau<-mimTau^2
diag(mimTau)<-0
maxi<-0.999999
mimTau[which(mimTau>maxi)]<-maxi
mimTau <--0.5*log(1-mimTau)
mimTau[mimTau<0]<-0
netTau <- mrnet(mimTau)

mimTau2Data <- readMat(sprintf("/Users/kiran/data/netbenchmark/matlab_outputs/syntren300_%d_tauoutput.mat",syntrenIdx))
mimTau2 = mimTau2Data$R
mimTau2 <- mimTau2^2
diag(mimTau2)<-0
maxi<-0.999999
mimTau2[which(mimTau2>maxi)]<-maxi
mimTau2 <--0.5*log(1-mimTau2)
mimTau2[mimTau2<0]<-0
netTau2 <- mrnet(mimTau2)

# make sure no N/A in the networks
if(sum(is.na(netCim)>0)){
  netCim[is.na(netCim)] <- 0
}
if(sum(is.na(netTau)>0)){
  netTau[is.na(netTau)] <- 0
}

ngenes <- dim(true.net)[1]
npos <- sum(true.net) #number of true links in the network
nlinks <- ngenes^2-ngenes #number of posible links in the network
no.topedges <- 20
if(sym){
  nlinks <- nlinks/2
}
no.edges <- round(nlinks*no.topedges/100)

rTau  <- evaluate(netTau,true.net,extend=no.edges,sym=TRUE)

### BEFORE ASSIGNING COL NAMES
# rTau2 <- evaluate(netTau2,true.net,extend=no.edges,sym=TRUE)
# rCim  <- evaluate(netCim,true.net,extend=no.edges,sym=TRUE)
# tauCalcLen = min(length(rTau[,"TP"]),length(rTau2[,"TP"]))
# print("***** BEFORE ASSIGNMENT *****")
# print(sprintf("(tau-tau2) = %d %d %d %d",
#       sum(rTau[1:tauCalcLen,"TP"]-rTau2[1:tauCalcLen,"TP"]),
#       sum(rTau[1:tauCalcLen,"FP"]-rTau2[1:tauCalcLen,"FP"]),
#       sum(rTau[1:tauCalcLen,"FP"]-rTau2[1:tauCalcLen,"TN"]),
#       sum(rTau[1:tauCalcLen,"FP"]-rTau2[1:tauCalcLen,"FN"])))
# print(sprintf("(tau-cim) = %d %d %d %d",
#       sum(rTau[1:tauCalcLen,"TP"]-rCim[1:tauCalcLen,"TP"]),
#       sum(rTau[1:tauCalcLen,"FP"]-rCim[1:tauCalcLen,"FP"]),
#       sum(rTau[1:tauCalcLen,"FP"]-rCim[1:tauCalcLen,"TN"]),
#       sum(rTau[1:tauCalcLen,"FP"]-rCim[1:tauCalcLen,"FN"])))

### AFTER ASSIGNING COL NAMES
colnames(netTau2) <- colnames(true.net)
rownames(netTau2) <- rownames(true.net)
colnames(netCim) <- colnames(true.net)
rownames(netCim) <- rownames(true.net)
rTau2 <- evaluate(netTau2,true.net,extend=no.edges,sym=TRUE)
rCim  <- evaluate(netCim,true.net,extend=no.edges,sym=TRUE)
tauCalcLen = min(length(rTau[,"TP"]),length(rTau2[,"TP"]))
cimCalcLen = min(length(rTau[,"TP"]),length(rCim[,"TP"]))
print("***** AFTER ASSIGNMENT *****")
print(sprintf("(tau-tau2) = %d %d %d %d",
      max(rTau[1:tauCalcLen,"TP"]-rTau2[1:tauCalcLen,"TP"]),
      max(rTau[1:tauCalcLen,"FP"]-rTau2[1:tauCalcLen,"FP"]),
      max(rTau[1:tauCalcLen,"TN"]-rTau2[1:tauCalcLen,"TN"]),
      max(rTau[1:tauCalcLen,"FN"]-rTau2[1:tauCalcLen,"FN"])))
print(sprintf("(tau-cim) = %d %d %d %d",
      max(rTau[1:cimCalcLen,"TP"]-rCim[1:cimCalcLen,"TP"]),
      max(rTau[1:cimCalcLen,"FP"]-rCim[1:cimCalcLen,"FP"]),
      max(rTau[1:cimCalcLen,"TN"]-rCim[1:cimCalcLen,"TN"]),
      max(rTau[1:cimCalcLen,"FN"]-rCim[1:cimCalcLen,"FN"])))


# par(mfrow=c(2,2))
# tauCalcLen = min(length(rTau[,"TP"]),length(rTau2[,"TP"]))
# plot(rTau[1:finalLen,"TP"]-rTau2[1:finalLen,"TP"],ylab="TP(tau-tau2)")
# plot(rTau[1:finalLen,"FP"]-rTau2[1:finalLen,"FP"],ylab="FP(tau-tau2)")
# plot(rTau[1:finalLen,"TN"]-rTau2[1:finalLen,"TN"],ylab="TN(tau-tau2)")
# plot(rTau[1:finalLen,"FN"]-rTau2[1:finalLen,"FN"],ylab="FN(tau-tau2)")

# par(mfrow=c(2,2))
# cimCalcLen = min(length(rTau[,"TP"]),length(rCim[,"TP"]))
# plot(rTau[1:finalLen,"TP"]-rCim[1:finalLen,"TP"],ylab="TP(tau-cim)")
# plot(rTau[1:finalLen,"FP"]-rCim[1:finalLen,"FP"],ylab="FP(tau-cim)")
# plot(rTau[1:finalLen,"TN"]-rCim[1:finalLen,"TN"],ylab="TN(tau-cim)")
# plot(rTau[1:finalLen,"FN"]-rCim[1:finalLen,"FN"],ylab="FN(tau-cim)")

# diffNet <- netCim-netTau
# which(diffNet!=0,arr.ind = T)

diffTau <- abs(mimTau-mimTau2)
which(diffTau>0.01,arr.ind=T)

diffNet <- abs(netTau-netTau2)
which(diffNet>0.01,arr.ind=T)

diffCimNet <- abs(netTau-netCim)
which(diffCimNet>0.1,arr.ind=T)

#image(mimCim-mimTau)
#library(corrplot)
#corrplot(mimTau, method="circle")
# plot(rTau[1:10000,"TP"],rCim[1:10000,"TP"])
# plot(rTau[1:10000,"FP"],rCim[1:10000,"FP"])
# plot(rTau[1:10000,"TN"],rCim[1:10000,"TN"])
# plot(rTau[1:10000,"FN"],rCim[1:10000,"FN"])

# tau_top20.aupr <- netbenchmark_custom.data(methods=c("Kendall"),data = ds,
#                                            true.net=true.net,plot=FALSE,verbose=FALSE)

# fullPathIn contains the path to the Matlab Syntren outputs
cim_top20.aupr <- netbenchmark_custom.data(methods=c("CIM"),data = fullPathIn,
                                           true.net=true.net,plot=FALSE,verbose=FALSE)
