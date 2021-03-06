function val = cm(X,Y,Z)

% val = cm(X,Y,Z,param) computes the CM statistic between two conditional CDFs
% Input: X, Y and Z are column vectors of length n, X is [nx1], number of columns can vary for Y and Z
% Output: val is the value of CM statistic
% Author: Sohan Seth (sohan@cnel.ufl.edu)

Kxx = grammat(X,'idfs');
Kyy = grammat(Y,'idfs');
Kzz = grammat(Z,'idfs');
Kyz = Kyy .* Kzz;
val = mean(mean(Kxx .* Kyz) .* mean(Kzz) - mean(Kxx .* Kzz) .* mean(Kyz)).^2;