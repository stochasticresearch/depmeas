## Combine the results from run_netbenchmark_tests.R for CIM and the remaining
rm(list = ls())
cat("\014")

##arrange df vars by position
##'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}

if(Sys.info()["sysname"]=="Darwin") {  # mac
  rResultsRepo_cnr = "/Users/kiran/data/netbenchmark/r_outputs"
  combinedResultsRepo_cnr = "/Users/kiran/data/netbenchmark/combined_outputs"
} else {  # unix
  rResultsRepo_cnr = "/home/kiran/data/netbenchmark/r_outputs"
  combinedResultsRepo_cnr = "/home/kiran/data/netbenchmark/combined_outputs"
}

numDatasetsToProcess = 150

# for(ds in c("syntren300","rogers1000","syntren1000","gnw1565","gnw2000")) {
for(ds in c("syntren300")) {
  # create an empty dataframe for this datasource, we merge each
  # individual result into this overall dataframe, and store the
  # merged output
  merged_df = data.frame()
  for(i in seq(numDatasetsToProcess)) {
    # load the R estimators data
    rf = file.path(rResultsRepo_cnr,sprintf("%s_%d.Rdata",ds,i))
    load(rf)
    # load the Matlab estimator data
    mf = file.path(rResultsRepo_cnr,sprintf("%s_%d_matlab.Rdata",ds,i))
    load(mf)
    # merge
    top20.aupr$`AUPRtop20%`$CIM   = cim_top20.aupr$`AUPRtop20%`$CIM
    top20.aupr$`AUPRtop20%`$KNN1  = knn1_top20.aupr$`AUPRtop20%`$MatlabMI
    top20.aupr$`AUPRtop20%`$KNN6  = knn6_top20.aupr$`AUPRtop20%`$MatlabMI
    top20.aupr$`AUPRtop20%`$KNN20 = knn20_top20.aupr$`AUPRtop20%`$MatlabMI
    top20.aupr$`AUPRtop20%`$vME   = vme_top20.aupr$`AUPRtop20%`$MatlabMI
    top20.aupr$`AUPRtop20%`$AP    = ap_top20.aupr$`AUPRtop20%`$MatlabMI
    
    merged_df <- rbind(merged_df, top20.aupr$`AUPRtop20%`)
  }
  # store off merged data
  merged_df <- arrange.vars(merged_df, c("CIM"=1,"Kendall"=2,"KNN1"=3,"KNN6"=4,"KNN20"=5,"vME"=6,"AP"=7))
  print(merged_df)
  fOut = file.path(combinedResultsRepo_cnr,sprintf("%s_combined.csv",ds))
  write.table(merged_df,fOut,sep=",",row.names=FALSE)
}
