#**************************************************************************
#* 
#* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
#*
#* This program is free software: you can redistribute it and/or modify
#* it under the terms of the GNU General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#*
#* You should have received a copy of the GNU General Public License
#* along with this program.  If not, see <http://www.gnu.org/licenses/>.
#**************************************************************************

rm(list = ls())
cat("\014")

# Figure out how to use the netbenchmark library
library(netbenchmark)
library(minet)
writeMat <- R.matlab::writeMat
readMat <- R.matlab::readMat

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

Pearson <- function(data){ 
  mim <- build.mim(data, estimator = "pearson")
  net <- mrnet(mim)
} 

CIM <- function(data) {
  fname <- paste(tempfile(), '.mat', sep='')
  writeMat(fname, dat=data)
  
  # matlab script location
  scriptLocation = '/home/kiran/stochasticresearch'
  scriptName = 'run_paircim_Rinterface.sh'
  MATLAB_runtime_path = '/usr/local/MATLAB/MATLAB_Runtime/v91'
  script = sprintf('%s %s %s', file.path(scriptLocation,scriptName), MATLAB_runtime_path, fname)
  
  # run the matlab script
  system(script)
  
  # read in the matlab data and convert to R format
  outputFname = sprintf('%s.matlab',fname)
  cimPairwise <- read.csv(outputFname,header=FALSE)
  colnames(cimPairwise) = colnames(data)
  rownames(cimPairwise) = colnames(data)
  
  # convert to matrix and build MRnet network
  mim <- data.matrix(cimPairwise)
  net <- mrnet(mim)
}

# we choose to run the datasources separately so that we can compare the performance for each dataset individually
# rather than an aggregate
comp <- netbenchmark(datasources.names="rogers1000", 
                     methods=c("CIM", "MIEmpirical", "MI.SG", "Kendall", "Pearson"),verbose=FALSE,
                     seed = 123) 
# save the results
fnameOut = "/home/kiran/ownCloud/PhD/sim_results/mrnet_rogers1000_cim.RData"
save(list = ls(all.names = TRUE), file = fnameOut, envir = .GlobalEnv)