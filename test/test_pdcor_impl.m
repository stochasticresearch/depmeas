%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
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

%% A script which tests the differences between the R and Matlab implementations
% of pdcorr

clear;
clc;

nsim = 100;
M = 500;
thresh = 0.001;
pval_thresh = 0.01;
numR = 2000;

pdcov_mse_vec = zeros(39, nsim);
pdcov_pval_mse_vec = zeros(39, nsim);
pdcor_mse_vec = zeros(39, nsim);
pdcor_pval_mse_vec = zeros(39, nsim);

% save the data
if(ispc)
    fname = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\test_pdcor_impl.txt';
elseif(ismac)
    fname = '/Users/Kiran/ownCloud/PhD/sim_results/independence/test_pdcor_impl.txt';
else
    fname = '/home/kiran/ownCloud/PhD/sim_results/independence/test_pdcor_impl.txt';
end
fID = fopen(fname, 'wt');

for ii=1:nsim
    fprintf(fID, '**************************** ii=%d ****************************\n', ii);
    
    % generate uncorrelated random variables
    x = rand(M,1); y = rand(M,1); z = rand(M,1);
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
        
    % generate correlated random variables
    xyz = mvnrnd([0 0 0], [1 0 .3; 0 1 -.1; .3 -.1 1], M);
    x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    xyz = mvnrnd([0 0 0], [1 0 .2; 0 1 -.8; .2 -.8 1], M);
    x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % generate ^ structure of various types
    
    % LINEAR ^-STRUCTURE
    x = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*x + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*x + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*x + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % QUADRATIC ^-STRUCTURE
    x = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    y = gamma*(x-0.5).^2 + (1-gamma)*eps;
    z = gamma*(x-0.5).^2 + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    y = gamma*(x-0.5).^2 + (1-gamma)*eps;
    z = gamma*(x-0.5).^2 + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    y = gamma*(x-0.5).^2 + (1-gamma)*eps;
    z = gamma*(x-0.5).^2 + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % SINUSOIDAL ^-STRUCTURE
    x = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    y = gamma*sin(4*pi*x) + (1-gamma)*eps;
    z = gamma*cos(4*pi*x) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    y = gamma*sin(4*pi*x) + (1-gamma)*eps;
    z = gamma*cos(4*pi*x) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    y = gamma*sin(4*pi*x) + (1-gamma)*eps;
    z = gamma*cos(4*pi*x) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % nthroot ^-STRUCTURE
    x = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    y = gamma*nthroot(x,4) + (1-gamma)*eps;
    z = gamma*nthroot(x,4) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    y = gamma*nthroot(x,4) + (1-gamma)*eps;
    z = gamma*nthroot(x,4) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    y = gamma*nthroot(x,4) + (1-gamma)*eps;
    z = gamma*nthroot(x,4) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % LINEAR/QUADRATIC ^-STRUCTURE
    x = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*(x-0.5).^2 + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*(x-0.5).^2 + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*(x-0.5).^2 + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % QUADRATIC/SINUSOIDAL ^-STRUCTURE
    x = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    y = gamma*(x-0.5).^2 + (1-gamma)*eps;
    z = gamm*cos(4*pi*x) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    y = gamma*(x-0.5).^2 + (1-gamma)*eps;
    z = gamm*cos(4*pi*x) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    x = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    y = gamma*(x-0.5).^2 + (1-gamma)*eps;
    z = gamm*cos(4*pi*x) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % generate v structure of various types
    
    % LINEAR V-Structure Testing
    y = rand(M,1); z = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    x = gamma*(y+z) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    x = gamma*(y+z) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    x = gamma*(y+z) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % QUADRATIC V-Structure Testing
    y = rand(M,1); z = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    x = gamma*((y-0.5).^2 + (z-0.5).^2) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    x = gamma*((y-0.5).^2 + (z-0.5).^2) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    x = gamma*((y-0.5).^2 + (z-0.5).^2) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % sinusoidal V-Structure Testing
    y = rand(M,1); z = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    x = gamma*(sin(4*pi*y)+cos(4*pi*z)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    x = gamma*(sin(4*pi*y)+cos(4*pi*z)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    x = gamma*(sin(4*pi*y)+cos(4*pi*z)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % nthroot V-Structure Testing
    y = rand(M,1); z = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    x = gamma*(nthroot(y,4)+nthroot(z,4)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    x = gamma*(nthroot(y,4)+nthroot(z,4)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    x = gamma*(nthroot(y,4)+nthroot(z,4)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % linear + quadratic V-Structure Testing
    y = rand(M,1); z = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    x = gamma*(y + (z-0.5).^2) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    x = gamma*(y + (z-0.5).^2) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    x = gamma*(y + (z-0.5).^2) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    % quadratic + sinusoidal V-Structure Testing
    y = rand(M,1); z = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    x = gamma*((y-0.5).^2 + cos(4*pi*z)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    x = gamma*((y-0.5).^2 + cos(4*pi*z)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    x = gamma*((y-0.5).^2 + cos(4*pi*z)) + (1-gamma)*eps;
    [pdcov_r, pdcov_pval_r] = pdcov_R(x,y,z, numR);
    [pdcor_r, pdcor_pval_r] = pdcorr_R(x,y,z, numR);
    [pdcov_matlab, pdcor_matlab, pdcov_pval_matlab, pdcor_pval_matlab] = pdcov(x,y,z, numR);
    fprintf(fID, 'pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf(fID, 'pdcov_pval_r=%0.06f pdcov_pval_matlab=%0.06f\n', pdcov_pval_r, pdcov_pval_matlab);
    fprintf(fID, 'pdcor_r=%0.06f pdcor_matlab=%0.06f\n', pdcor_r, pdcor_matlab);
    fprintf(fID, 'pdcor_pval_r=%0.06f pdcor_pval_matlab=%0.06f\n', pdcor_pval_r, pdcor_pval_matlab);
    fprintf(fID, '\n');
    
    storeIdx = storeIdx + 1;
    pdcov_mse_vec(storeIdx,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(storeIdx,ii) = (pdcov_pval_r-pdcov_pval_matlab).^2;
    pdcor_mse_vec(storeIdx,ii) = (pdcor_r-pdcor_matlab).^2;
    pdcor_pval_mse_vec(storeIdx,ii) = (pdcor_pval_r-pdcor_pval_matlab).^2;
    
end

pdcov_mse = mean(pdcov_mse_vec,2);
pdcov_pval_mse = mean(pdcov_pval_mse_vec,2);
pdcor_mse = mean(pdcor_mse_vec,2);
pdcor_pval_mse = mean(pdcor_pval_mse_vec,2);

fprintf(fID, '>>>>>>>>>>>>>>>>>>>>>>>> RESULTS <<<<<<<<<<<<<<<<<<<<<\n');
fprintf(fID, 'pdcov_mse=%0.03f\n', pdcov_mse);
fprintf(fID, 'pdcov_pval_mse=%0.03f\n', pdcov_pval_mse);
fprintf(fID, 'pdcor_mse=%0.03f\n', pdcor_mse);
fprintf(fID, 'pdcor_pval_mse=%0.03f\n', pdcor_pval_mse);
fprintf(fID, '>>>>>>>>>>>>>>>>>>>>>>>>*********<<<<<<<<<<<<<<<<<<<<<\n');

fclose(fID);
