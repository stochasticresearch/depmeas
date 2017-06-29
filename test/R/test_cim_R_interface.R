# test out the matlab/r interface

rm(list = ls())
cat("\014")

library(minet)
writeMat <- R.matlab::writeMat
readMat <- R.matlab::readMat

data(syn.data)
data_matrix = data.matrix(syn.data)
fname <- paste(tempfile(), '.mat', sep='')
writeMat(fname, dat=data_matrix)

# matlab script location
scriptLocation = '/home/kiran/stochasticresearch'
scriptName = 'run_paircim_Rinterface.sh'
MATLAB_runtime_path = '/usr/local/MATLAB/MATLAB_Runtime/v91'
script = sprintf('%s %s %s', file.path(scriptLocation,scriptName), MATLAB_runtime_path, fname)

system(script)

# read in the matlab data and convert to R format
outputFname = sprintf('%s.matlab',fname)
cimPairwise <- read.csv(outputFname,header=FALSE)
colnames(cimPairwise) = colnames(data_matrix)
# convert from dataframe to matrix
cim_mim <- data.matrix(cimPairwise)
cim_net <- mrnet(cim_mim)

mim <- build.mim(data_matrix, estimator = "mi.empirical", disc = "equalfreq")
net <- mrnet(mim)