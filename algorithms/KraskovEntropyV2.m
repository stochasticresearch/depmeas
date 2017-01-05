function [ H ] = KraskovEntropyV2( X, k, varargin )
%KraskovEntropy computes the Kraskov estimator of the Shannon entropy of
%the dataset X, using the k-nearest neighbours sample. By default it uses
%the euclidean distance

%   X: dataset (num samples x num variables)
%   k: number of nearest neighbours
%   varargin: distance (default 'euclidean')

% Reference:
% 1) Kraskov, A.; Stogbauer, H.; Grassberger, P. "Estimating mutual information".
%    Phys. Rev. E Stat. Nonlin. Soft. Matter Phys. 2004, 69, 1–15.
%
% 2) Veselkov K. A. et al., "A Metabolic Entropy Approach for Measurements
%    of Systemic Metabolic Disruptions in Patho-Physiological States".
%    Journal of proteome research 9 (7), 3537-3544
%    doi:10.1021/pr1000576

nSamples = size(X, 1);
nVariables = size(X, 2);

if k > nSamples
    error('k must be smaller than number of samples.');
end

if nnz(isnan(X)) > 0
    error('X cannot contain NaN values.');
end

if nargin == 2
    distance = 'euclidean';
else
    distance = varargin{1};
end

[~, distEpsilon] = knnsearch(X, X, 'K', k+1, 'distance', distance);
distEpsilon(:, 1) = []; % remove self-distance (all zeros)

Cd = pi^(nVariables/2) / gamma(1 + nVariables/2) / 2^nVariables;
H = psi(nSamples) - psi(k) + log(Cd) + (nVariables / nSamples) * sum(log(2*distEpsilon(:,k)));

end

