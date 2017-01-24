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

Spearman <- function(data){ 
  mim <- build.mim(data, estimator = "spearman", disc = "none")
  net <- mrnet(mim)
} 

Kendall <- function(data){ 
  mim <- build.mim(data, estimator = "kendall", disc = "none")
  net <- mrnet(mim)
} 

MI <- function(data){ 
  mim <- build.mim(data, estimator = "mi.empirical", disc = "equalfreq")
  net <- mrnet(mim)
} 


comp <- netbenchmark(datasources.names="all", 
                     methods=c("Kendall","MI"),verbose=FALSE,
                     seed = 123) 
# save the results
save(list = ls(all.names = TRUE), file = "/home/kiran/ownCloud/PhD/sim_results/mrnet_all.RData", envir = .GlobalEnv)
aupr <- comp[[1]][,-(1:2)]
#make the name look prety 
library("tools") 
colnames(aupr) <- sapply(colnames(aupr),file_path_sans_ext) 
boxplot(aupr, main="All Data", ylab=expression('AUPR'[20]))