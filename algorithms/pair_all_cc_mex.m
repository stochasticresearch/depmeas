function [R_cim,R_mi_knn1,R_mi_knn6,R_mi_knn20,R_mi_vme,R_mi_ap] = ...
    pair_all_cc_mex( X )
%PAIRCIM - computes pairwise dependency metrics of a given vector
% Inputs:
%  X - a matrix of observations from which pairwise CIM statistics are
%  computed.
% Outputs:
%  R - a matrix which contains pairwise correlations of all columns in X
%**************************************************************************
%*                                                                        *
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

n = size(X,2);      % the dimensionality of the dataset
R_cim = zeros(n,n);
R_mi_knn1  = zeros(n,n);
R_mi_knn6  = zeros(n,n);
R_mi_knn20 = zeros(n,n);
R_mi_vme   = zeros(n,n);
R_mi_ap    = zeros(n,n);

minScanIncr = 0.015625;

for ii=1:n
    x = X(:,ii);
    parfor (jj=ii+1:n, getParforArg())
        y = X(:,jj);        
        
        cimVal = cim_v8a_cc_mex(x,y,minScanIncr);
        knn1Val = KraskovMI_cc_mex(x,y,1);
        knn6Val = KraskovMI_cc_mex(x,y,6);
        knn20Val = KraskovMI_cc_mex(x,y,20);
        vmeVal = vmeMI_interface(x,y);
        apVal  = apMI_interface(x,y);
        
        R_cim(ii,jj)      = cimVal;
        R_mi_knn1(ii,jj)  = knn1Val;
        R_mi_knn6(ii,jj)  = knn6Val;
        R_mi_knn20(ii,jj) = knn20Val;
        R_mi_vme(ii,jj)   = vmeVal;
        R_mi_ap(ii,jj)    = apVal;
    end
end

% assign the lower triangle by copying the upper-triangle of the R matrix
% make matrix symmetric.  The below works because R is initialized to
% zeros.
R_cim=R_cim+R_cim';
R_cim(1:n+1:n*n) = 0;   % set the diagonal to 0, to match R's MINET package

R_mi_knn1=R_mi_knn1+R_mi_knn1';
R_mi_knn1(1:n+1:n*n) = 0;   % set the diagonal to 0, to match R's MINET package

R_mi_knn6=R_mi_knn6+R_mi_knn6';
R_mi_knn6(1:n+1:n*n) = 0;   % set the diagonal to 0, to match R's MINET package

R_mi_knn20=R_mi_knn20+R_mi_knn20';
R_mi_knn20(1:n+1:n*n) = 0;   % set the diagonal to 0, to match R's MINET package

R_mi_vme=R_mi_vme+R_mi_vme';
R_mi_vme(1:n+1:n*n) = 0;   % set the diagonal to 0, to match R's MINET package

R_mi_ap=R_mi_ap+R_mi_ap';
R_mi_ap(1:n+1:n*n) = 0;   % set the diagonal to 0, to match R's MINET package

end