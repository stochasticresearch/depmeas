rm(list = ls())
cat("\014")

library(R.matlab)
library(distr)
library(distrEx)

betaDistInfo <- readMat('/home/kiran/ownCloud/PhD/sim_results/independence/cim_nullDistribution_beta.mat')

ii = 1
hd_vec <- vector(mode="numeric", length=length(betaDistInfo$M.vec))
for (M in betaDistInfo$M.vec) {
  alpha1 = 0.2;
  alpha2 = 1-alpha1;
  M1 = trunc(alpha1*M);
  M2 = trunc(alpha1*M);
  d1 = Norm(mean=0,sd=2*(2*M1+5)/(9*M1*(M1-1)));
  d2 = Norm(mean=0,sd=2*(2*M2+5)/(9*M2*(M2-1)));
  trueNullDist_approx = alpha1*abs(d1)+alpha1*abs(d2);
  
  empiricalNullDist_approx = Beta(shape1 = betaDistInfo$alphaVecContinuous[ii], shape2 = betaDistInfo$betaVecContinuous[ii])

  hd_vec[ii] = HellingerDist(trueNullDist_approx,empiricalNullDist_approx)
  
  ii = ii +1;
}
plot(betaDistInfo$M.vec,hd_vec,xlab="M (Sample Size)", ylab="Hellinger Distance",log="x")
writeMat('/home/kiran/ownCloud/PhD/sim_results/independence/cim_nullDistribution_hd.mat', M_vec = betaDistInfo$M.vec, hd_vec=hd_vec)