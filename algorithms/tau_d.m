function [ m ] = tau_d( X )
%TAU_D Estimates multivariate Kendall's tau - as defined by Joe (1990) in
%the paper: Multivariate Concordance (1990)

U = pobs(X);
% a vectorized but naive way to compute! -- research
% better ways perhaps ...
sumVal = 0;
for ii=1:M
    sumVal = sumVal + sum(sum(U(ii,:)<U(ii+1:end,:),2)==d);
end
m = (2^d/(M*(M-1))*sumVal-1)*(1/(2^(d-1)-1));


end

