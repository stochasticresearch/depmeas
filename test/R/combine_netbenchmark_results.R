## Combine the results from run_netbenchmark_tests.R for CIM and the remaining
rm(list = ls())
cat("\014")

if(Sys.info()["sysname"]=="Darwin") {  # mac
  rResultsRepo_cnr = "/Users/kiran/data/netbenchmark/r_outputs"
  combinedResultsRepo_cnr = "/Users/kiran/data/netbenchmark/combined_outputs"
} else {  # unix
  rResultsRepo_cnr = "/home/kiran/data/netbenchmark/r_outputs"
  combinedResultsRepo_cnr = "/home/kiran/data/netbenchmark/combined_outputs"
}

# numDatasetsToProcess = 150
numDatasetsToProcess = 3

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
    mf = file.path(rResultsRepo_cnr,sprintf("%s_%d_CIM.Rdata",ds,i))
    load(mf)
    # merge
    top20.aupr$`AUPRtop20%`$CIM = cim_top20.aupr$`AUPRtop20%`$CIM
    merged_df <- rbind(merged_df, top20.aupr$`AUPRtop20%`)
  }
  # store off merged data
  print(merged_df)
  fOut = file.path(combinedResultsRepo_cnr,sprintf("%s_combined.csv",ds))
  write.table(merged_df,fOut,sep=",",row.names=FALSE)
}
