## Convert the Rdata files to Matlab files, so we can get the MIM for the CIM metric

rm(list = ls())
cat("\014")

writeMat <- R.matlab::writeMat
readMat <- R.matlab::readMat

availableDataSources = c("rogers1000","syntren300","syntren1000","gnw1565","gnw2000")
dataRepo = "/home/kiran/data/netbenchmark/inputs/"

for(ds in availableDataSources) {
  inputFname = sprintf("%s.Rdata", ds)
  fullPathIn = file.path(dataRepo,inputFname)
  load(file=fullPathIn)
  
  for(i in seq_along(data.list)){
    outputFname = sprintf("%s_%d.mat",ds,i)
    fullPathOut = file.path(dataRepo,outputFname)
    
    sd <- data.list[[i]]
    writeMat(fullPathOut,data=sd)
  }
}