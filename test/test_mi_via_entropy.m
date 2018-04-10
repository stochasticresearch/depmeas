%% Estimate the theoretical mutual information between RV's of different dependencies using Copulas

delta_uv = 100;
u = linspace(0,1,delta_uv);
[u1,u2] = meshgrid(u,u);
    
tauVec = linspace(0.01,0.99,25);
MI_estimate = zeros(4,length(tauVec));
for tauIdx=1:length(tauVec)
    tau = tauVec(tauIdx);
    % for each copula type, we compute the mutual information of *any* two
    % marginal random variables 
    theta = copulaparam('Gaussian',tau,'type','kendall');
    c_gauss = copulapdf('Gaussian',[u1(:),u2(:)],theta);
    cc_gauss = reshape(c_gauss,delta_uv,delta_uv);
    
    theta = copulaparam('Frank',tau,'type','kendall');
    c_frank = copulapdf('Frank',[u1(:),u2(:)],theta);
    cc_frank = reshape(c_frank,delta_uv,delta_uv);
    
    theta = copulaparam('Gumbel',tau,'type','kendall');
    c_gumbel = copulapdf('Gumbel',[u1(:),u2(:)],theta);
    cc_gumbel = reshape(c_gumbel,delta_uv,delta_uv);
    
    theta = copulaparam('Clayton',tau,'type','kendall');
    c_clayton = copulapdf('Clayton',[u1(:),u2(:)],theta);
    cc_clayton = reshape(c_clayton,delta_uv,delta_uv);
    
    
    % compute the MI for each copula, estimated based on entropy
    % computation
    MI_estimate(1,tauIdx) = -1*entropy2d(cc_gauss,u1,u2);
    MI_estimate(2,tauIdx) = -1*entropy2d(cc_frank,u1,u2);
    MI_estimate(3,tauIdx) = -1*entropy2d(cc_gumbel,u1,u2);
    MI_estimate(4,tauIdx) = -1*entropy2d(cc_clayton,u1,u2);
    
end

%% Compute MI between Hybrid X&Y using the conditional entropy formula
% MI(X,Y) = H(Y) - H(Y|X), where Y is discrete, and X can be either
% discrete or continuous

clear;
clc;

leftSkewData = pearsrnd(0,1,-1,3,5000,1);
rightSkewData = pearsrnd(0,1,1,3,5000,1);

[fLeftSkew,xiLeftSkew] = emppdf(leftSkewData,0);
FLeftSkew = empcdf(xiLeftSkew,0);
leftSkewContinuousDistInfo = rvEmpiricalInfo(xiLeftSkew,fLeftSkew,FLeftSkew,0);
[fRightSkew,xiRightSkew] = emppdf(rightSkewData,0);
FRightSkew = empcdf(xiRightSkew,0);
rightSkewContinuousDistInfo = rvEmpiricalInfo(xiRightSkew,fRightSkew,FRightSkew,0);
distObj = makedist('Normal');

cop = 'Gaussian';
tau = 0.6;
M = 500;
iTau = copulaparam(cop,tau);

U = copularnd(cop,iTau,M);
X = icdf(distObj,U(:,1));
distObj = makedist('Multinomial','probabilities',[0.5,0.5]);
Y = icdf(distObj,U(:,2));
[X,Y] = pobs_sorted_cc(X,Y);

I_XY = discrete_entropy(Y) - conditional_entropy(X,Y,1);