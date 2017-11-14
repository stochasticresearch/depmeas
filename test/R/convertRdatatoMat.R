## Convert the Rdata files to Matlab files, so we can get the MIM for the CIM metric

rm(list = ls())
cat("\014")

writeMat <- R.matlab::writeMat
readMat <- R.matlab::readMat

#availableDataSources = c("syntren300","rogers1000","syntren1000","gnw1565","gnw2000")
availableDataSources = c("syntren300")
localNoiseVec = c(0,10,20,30,40,50)
globalNoiseVec = c(0,10,20,30,40,50)
dataRepo = "/data/netbenchmark/inputs/"

for(ds in availableDataSources) {
    for(gn in globalNoiseVec) {
        for(ln in localNoiseVec) {
            inputFname = sprintf("%s.Rdata", ds)
            subFolder = sprintf('gn_%d_ln_%d',gn,ln) 
            fullPathIn = file.path(dataRepo,subFolder,inputFname)
            load(file=fullPathIn)
            
            for(i in seq_along(data.list)){
                outputFname = sprintf("%s_%d.mat",ds,i)
                dir.create(file.path(dataRepo,subFolder), showWarnings = FALSE)
                fullPathOut = file.path(dataRepo,subFolder,outputFname)
                
                sd <- data.list[[i]]
                writeMat(fullPathOut,data=sd)
            }
        }
    }
}