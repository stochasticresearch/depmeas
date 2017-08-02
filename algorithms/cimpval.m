function [pval] = cimpval(cimVal, M, varargin)
%CIMPVAL - computes the p-value of a given CIM statistic and the number of
%samples upon which that CIM statistic was calculated against independence
%hypothesis
% Inputs:
%  cimVal - the CIM measure
%  M - the sample size used to compute this CIM measure
%  varargin{1} - type of data for which this CIM metric was computed,
%                continuous, discrete, hybrid1, or hybrid2.  continuous is
%                default
% values are hard-coded, look at the script cim_runPower which has a
% execution-cell inside it that generates these vectors!

% loads the vectors: alphaVecContinuous  betaVecContinuous,
%                    alphaVecHybrid1     betaVecHybrid1
%                    alphaVecHybrid2     betaVecHybrid2
%                    alphaVecDiscrete    betaVecDiscrete
load('cimpval_data.mat');

if(length(varargin)>=1)
    selector = varargin{1};
else
    selector = 'continuous';
end

if(strcmpi(selector, 'continuous'))
    alphaVec = alphaVecContinuous;
    betaVec = betaVecContinuous;
elseif(strcmpi(selector, 'hybrid1'))
    alphaVec = alphaVecHybrid1;
    betaVec = betaVecHybrid1;
elseif(strcmpi(selector, 'hybrid2'))
    alphaVec = alphaVecHybrid2;
    betaVec = betaVecHybrid2;
else    % assume discrete
    alphaVec = alphaVecDiscrete;
    betaVec = betaVecDiscrete;
end

M_vec = [100:100:2500 5000 10000];

% iterpolate for each value of k, mu, sigma for the given M
alpha = interp1(M_vec, alphaVec, M, 'spline');
beta = interp1(M_vec, betaVec, M, 'spline');

pval = 1-betacdf(cimVal, alpha, beta);

end