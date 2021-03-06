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

%% Conditional Dependence Parametric test
clear;
clc;
dbstop if error;

rng(21);

significance_thresh = 0.05;

M = 500;
nsim = 500;

gammaVec = 0:0.1:1;

numDepTypes = 6;
cimResultsMat = zeros(length(gammaVec), numDepTypes, nsim);        % for comparision of acceptance rates
ccimResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
cmaSurrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);     % for comparision of aceptance rate
cmaResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pdcorrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pdcorrPValMat = zeros(length(gammaVec), numDepTypes, nsim);
pcorrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pcorrPValMat = zeros(length(gammaVec), numDepTypes, nsim);

cimResultsVec = zeros(1,nsim);
ccimResultsVec = zeros(1,nsim);
cmaSurrResultsVec = zeros(1,nsim);
cmaResultsVec = zeros(1,nsim);
pdcorrResultsVec = zeros(1,nsim);
pdcorrPValVec = zeros(1,nsim);
pcorrResultsVec = zeros(1,nsim);
pcorrPValVec = zeros(1,nsim);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for gammaIdx=1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    dispstat(sprintf('Computing for gamma=%0.02f',gamma),'keepthis', 'timestamp');
    for jj=1:numDepTypes
        dispstat(sprintf('DepType=%d',jj),'keepthis', 'timestamp');
        parfor ii=1:nsim
            %dispstat(sprintf('Simulating -- %0.02f %%', ii/nsim*100),'timestamp');
            Y = rand(M,1);
            Z = rand(M,1);
            eps = randn(M,1);
            
            switch(jj)
                case 1
                    X = gamma*(Y+Z) + (1-gamma)*eps;
                case 2
                    X = gamma*((Y-0.5).^2 + (Z-0.5).^2) + (1-gamma)*eps;
                case 3
                    X = gamma*(sin(4*pi*Y)+cos(4*pi*Z)) + (1-gamma)*eps;
                case 4
                    X = gamma*(nthroot(Y,4)+nthroot(Z,4)) + (1-gamma)*eps;
                case 5
                    X = gamma*(Y + (Z-0.5).^2) + (1-gamma)*eps;
                case 6
                    X = gamma*((Y-0.5).^2 + cos(4*pi*Z)) + (1-gamma)*eps;
            end
            
            cimVal = cim(Y,Z);        % for calculating the CI/CD test statistic
            ccimVal = ccim(Y,Z,X);
            [~, pdcorVal, ~, pdcor_pval] = pdcov(Y,Z,X);
            [pcorrVal,pcorr_pval] = partialcorr(Y,Z,X);

            data = struct();
            data.X = Y; data.Y = Z; data.Z = X;
            cmaVal = cassor(data);
            cmaSurrData = gensurr(data);
            surrData = struct();
            surrData.X = cmaSurrData.Y; surrData.Y = cmaSurrData.Z; surrData.Z = cmaSurrData.X;
            cmaSurrVal = cassor(surrData);
            
            cimResultsVec(ii) = cimVal;
            ccimResultsVec(ii) = ccimVal;
            
            cmaSurrResultsVec(ii) = cmaSurrVal;
            cmaResultsVec(ii) = cmaVal;
            
            pdcorrResultsVec(ii) = pdcorVal;
            pdcorrPValVec(ii) = pdcor_pval;
            
            pcorrResultsVec(ii) = pcorrVal;
            pcorrPValVec(ii) = pcorr_pval;
        end
        
        cimResultsMat(gammaIdx, jj, :) = cimResultsVec;
        ccimResultsMat(gammaIdx, jj, :) = ccimResultsVec;
        
        cmaSurrResultsMat(gammaIdx, jj, :) = cmaSurrResultsVec;
        cmaResultsMat(gammaIdx, jj, :) = cmaResultsVec;
        
        pdcorrResultsMat(gammaIdx, jj, :) = pdcorrResultsVec;
        pdcorrPValMat(gammaIdx, jj, :) = pdcorrPValVec;
        
        pcorrResultsMat(gammaIdx, jj, :) = pcorrResultsVec;
        pcorrPValMat(gammaIdx, jj, :) = pcorrPValVec;
    end
    % save the results in case the computer crashes :(
    if(ispc)
        save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\ccim_CD.mat');
    elseif(ismac)
        save('/Users/Kiran/ownCloud/PhD/sim_results/independence/ccim_CD.mat');
    elseif(isunix)
        save('/home/kiran/ownCloud/PhD/sim_results/independence/ccim_CD.mat');
    end

end

% To understand the acceptance for CMA, look at demo.m in the cmi folder to
% see how to use the "surrogate" versus the actual test statistic to see if
% there is an accept or reject.

% save the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\ccim_CD.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/ccim_CD.mat');
elseif(isunix)
    save('/home/kiran/ownCloud/PhD/sim_results/independence/ccim_CD.mat');
end

% compute acceptance rates
pdcorr_acceptance = pdcorrPValMat<significance_thresh;
pcorr_acceptance = pcorrPValMat<significance_thresh;
ccim_acceptance = ccimResultsMat>cimResultsMat;
cma_acceptance = cmaResultsMat>cmaSurrResultsMat;

% plot the dependence metric results versus gamma for each dep type
figure;

h1 = subplot(3,2,1);
depIdx = 1;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, meanccimAcceptanceRate, 'o-.', ...
     gammaVec, meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('Acceptance Rate', 'FontSize', 20);
title('(a)', 'FontSize', 20);
grid on;
h1.FontSize = 20; 

h2 = subplot(3,2,2);
depIdx = 2;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, meanccimAcceptanceRate, 'o-.', ...
     gammaVec, meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('Acceptance Rate', 'FontSize', 20);
title('(b)', 'FontSize', 20);
grid on;
h2.FontSize = 20;

h3 = subplot(3,2,3);
depIdx = 3;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, meanccimAcceptanceRate, 'o-.', ...
     gammaVec, meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('Acceptance Rate', 'FontSize', 20); 
title('(c)', 'FontSize', 20);
grid on;
h3.FontSize = 20;

h4 = subplot(3,2,4);
depIdx = 4;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, meanccimAcceptanceRate, 'o-.', ...
     gammaVec, meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('Acceptance Rate', 'FontSize', 20); 
title('(d)', 'FontSize', 20);
grid on;
h4.FontSize = 20;

h5 = subplot(3,2,5);
depIdx = 5;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, meanccimAcceptanceRate, 'o-.', ...
     gammaVec, meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('Acceptance Rate', 'FontSize', 20);
title('(e)', 'FontSize', 20);
grid on;
h5.FontSize = 20; 

h6 = subplot(3,2,6);
depIdx = 6;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, meanccimAcceptanceRate, 'o-.', ...
     gammaVec, meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('Acceptance Rate', 'FontSize', 20);
title('(f)', 'FontSize', 20);
grid on;
h6.FontSize = 20; 

legend('CCIM', 'PDCORR', 'PCORR', 'CMA');  % manually move this using the mouse to a
                                            % good location

%% Conditional Independence Test

clear;
clc;
dbstop if error;

rng(22);

M = 500;
nsim = 500;

significance_thresh = 0.05;
gammaVec = 0:0.1:1;

numDepTypes = 6;
cimResultsMat = zeros(length(gammaVec), numDepTypes, nsim);        % for comparision of acceptance rates
ccimResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
cmaSurrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);     % for comparision of aceptance rate
cmaResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pdcorrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pdcorrPValMat = zeros(length(gammaVec), numDepTypes, nsim);
pcorrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pcorrPValMat = zeros(length(gammaVec), numDepTypes, nsim);

cimResultsVec = zeros(1,nsim);
ccimResultsVec = zeros(1,nsim);
cmaSurrResultsVec = zeros(1,nsim);
cmaResultsVec = zeros(1,nsim);
pdcorrResultsVec = zeros(1,nsim);
pdcorrPValVec = zeros(1,nsim);
pcorrResultsVec = zeros(1,nsim);
pcorrPValVec = zeros(1,nsim);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for gammaIdx=1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    dispstat(sprintf('Computing for gamma=%0.02f',gamma),'keepthis', 'timestamp');
    for jj=1:numDepTypes
        dispstat(sprintf('DepType=%d',jj),'keepthis', 'timestamp');
        parfor ii=1:nsim
            %dispstat(sprintf('Simulating -- %0.02f %%', ii/nsim*100),'timestamp');
            X = rand(M,1);
            eps = randn(M,1);
            
            switch(jj)
                case 1
                    Y = gamma*X + (1-gamma)*eps;
                    Z = gamma*X + (1-gamma)*eps;
                case 2
                    Y = gamma*(X-0.5).^2 + (1-gamma)*eps;
                    Z = gamma*(X-0.5).^2 + (1-gamma)*eps;
                case 3
                    Y = gamma*sin(4*pi*X) + (1-gamma)*eps;
                    Z = gamma*cos(4*pi*X) + (1-gamma)*eps;
                case 4
                    Y = gamma*nthroot(X,4) + (1-gamma)*eps;
                    Z = gamma*nthroot(X,4) + (1-gamma)*eps;
                case 5
                    Y = gamma*X + (1-gamma)*eps;
                    Z = gamma*(X-0.5).^2 + (1-gamma)*eps;
                case 6
                    Y = gamma*(X-0.5).^2 + (1-gamma)*eps;
                    Z = gamma*cos(4*pi*X) + (1-gamma)*eps;
            end
            
            cimVal = cim(Y,Z);        % for calculating the CI/CD test statistic
            ccimVal = ccim(Y,Z,X);
            [~, pdcorVal, ~, pdcor_pval] = pdcov(Y,Z,X);
            [pcorrVal,pcorr_pval] = partialcorr(Y,Z,X);
            
            data = struct();
            data.X = Y; data.Y = Z; data.Z = X;
            cmaVal = cassor(data);
            cmaSurrData = gensurr(data);
            surrData = struct();
            surrData.X = cmaSurrData.Y; surrData.Y = cmaSurrData.Z; surrData.Z = cmaSurrData.X;
            cmaSurrVal = cassor(surrData);
            
            cimResultsVec(ii) = cimVal;
            ccimResultsVec(ii) = ccimVal;
            
            cmaSurrResultsVec(ii) = cmaSurrVal;
            cmaResultsVec(ii) = cmaVal;
            
            pdcorrResultsVec(ii) = pdcorVal;
            pdcorrPValVec(ii) = pdcor_pval;
            
            pcorrResultsVec(ii) = pcorrVal;
            pcorrPValVec(ii) = pcorr_pval;
        end
        
        cimResultsMat(gammaIdx, jj, :) = cimResultsVec;
        ccimResultsMat(gammaIdx, jj, :) = ccimResultsVec;
        
        cmaSurrResultsMat(gammaIdx, jj, :) = cmaSurrResultsVec;
        cmaResultsMat(gammaIdx, jj, :) = cmaResultsVec;
        
        pdcorrResultsMat(gammaIdx, jj, :) = pdcorrResultsVec;
        pdcorrPValMat(gammaIdx, jj, :) = pdcorrPValVec;
        
        pcorrResultsMat(gammaIdx, jj, :) = pcorrResultsVec;
        pcorrPValMat(gammaIdx, jj, :) = pcorrPValVec;
        
    end
    % save the results in case the computer crashes
    if(ispc)
        save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\ccim_CI.mat');
    elseif(ismac)
        save('/Users/Kiran/ownCloud/PhD/sim_results/independence/ccim_CI.mat');
    elseif(isunix)
        save('/home/kiran/ownCloud/PhD/sim_results/independence/ccim_CI.mat');
    end
end

% save the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\ccim_CI.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/ccim_CI.mat');
elseif(isunix)
    save('/home/kiran/ownCloud/PhD/sim_results/independence/ccim_CI.mat');
end

% compute acceptance rates
pdcorr_acceptance = pdcorrPValMat<significance_thresh;
pcorr_acceptance = pcorrPValMat<significance_thresh;
ccim_acceptance = ccimResultsMat>cimResultsMat;
cma_acceptance = cmaResultsMat>cmaSurrResultsMat;

figure;

h1 = subplot(3,2,1);
depIdx = 1;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, 1-meanccimAcceptanceRate, 'o-.', ...
     gammaVec, 1-meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, 1-meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, 1-meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('1-Acceptance Rate', 'FontSize', 20);
title('(a)', 'FontSize', 20);
grid on;
h1.FontSize = 20; 

h2 = subplot(3,2,2);
depIdx = 2;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, 1-meanccimAcceptanceRate, 'o-.', ...
     gammaVec, 1-meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, 1-meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, 1-meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('1-Acceptance Rate', 'FontSize', 20);
title('(b)', 'FontSize', 20);
grid on;
h2.FontSize = 20;

h3 = subplot(3,2,3);
depIdx = 3;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, 1-meanccimAcceptanceRate, 'o-.', ...
     gammaVec, 1-meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, 1-meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, 1-meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('1-Acceptance Rate', 'FontSize', 20); 
title('(c)', 'FontSize', 20);
grid on;
h3.FontSize = 20;

h4 = subplot(3,2,4);
depIdx = 4;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, 1-meanccimAcceptanceRate, 'o-.', ...
     gammaVec, 1-meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, 1-meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, 1-meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('1-Acceptance Rate', 'FontSize', 20); 
title('(d)', 'FontSize', 20);
grid on;
h4.FontSize = 20;

h5 = subplot(3,2,5);
depIdx = 5;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, 1-meanccimAcceptanceRate, 'o-.', ...
     gammaVec, 1-meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, 1-meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, 1-meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('1-Acceptance Rate', 'FontSize', 20);
title('(e)', 'FontSize', 20);
grid on;
h5.FontSize = 20; 

h6 = subplot(3,2,6);
depIdx = 6;
meanccimRes = mean(squeeze(ccimResultsMat(:,depIdx,:)),2);
meanccimAcceptanceRate = mean(squeeze(ccim_acceptance(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPdcorrAcceptanceRate = mean(squeeze(pdcorr_acceptance(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanPcorrAcceptanceRate = mean(squeeze(pcorr_acceptance(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
meanCmaAcceptanceRate = mean(squeeze(cma_acceptance(:,depIdx,:)),2);
plot(gammaVec, 1-meanccimAcceptanceRate, 'o-.', ...
     gammaVec, 1-meanPdcorrAcceptanceRate, 'x-.', ...
     gammaVec, 1-meanPcorrAcceptanceRate, 'd-.', ...
     gammaVec, 1-meanCmaAcceptanceRate, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('1-Acceptance Rate', 'FontSize', 20);
title('(f)', 'FontSize', 20);
grid on;
h6.FontSize = 20; 

legend('CCIM', 'PDCORR', 'PCORR', 'CMA');  % manually move this using the mouse to a
                                            % good location

%% Do a conditionally independent test

clear;
clc;

rng(123);

M = 500;
noise = 0.01;
alpha = 0.05;

% Generate data from     Y<--X-->Z
x = rand(M,1);
% y = x + noise*randn(M,1); z = x + noise*randn(M,1);
% y = 4*(x-0.5).^2 + noise*randn(M,1); z = 4*(x-0.5).^2 + noise*randn(M,1);
% y = sin(4*pi*x) + noise*randn(M,1); z = cos(4*pi*x) + noise*randn(M,1);
% y = nthroot(x, 4) + noise*randn(M,1); z = nthroot(x, 4) + noise*randn(M,1);
% y = x + noise*randn(M,1); z = 4*(x-0.5).^2 + noise*randn(M,1);
y = 4*(x-0.5).^2 + noise*randn(M,1); z = cos(4*pi*x) + noise*randn(M,1);

cim1 = cim(x, y);
cim2 = cim(x, z);
cim3 = cim(y, z);
% ccim conditions X on Y and Z, and sees how related Y and Z
% are to each other ... in this graphical model, they should be UNRELATED
% after the effect of X is removed... i.e. close to independent.  To see
% why, look at the graphical model, Y indep Z | X according to
% d-separation.  So if we condition upon X (i.e. remove teh effect of X on
% Y and Z separately), then we should get independence.
[ccimVal, RxAligned, RyAligned] = ccim(y,z,x);
[~, pdcorVal] = pdcov(y, z, x); pdcorVal = abs(pdcorVal);
partialCorrVal = abs(partialcorr(y,z,x));
data.X = y; data.Y = z; data.Z = x; cassorVal = abs(cassor(data));
hdVal = abs(hd(y, z, x));

fontSize = 20;

h1 = subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', fontSize); ylabel('y', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim1), 'FontSize', fontSize);
h1.FontSize = fontSize;

h2 = subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim3), 'FontSize', fontSize);
h2.FontSize = fontSize;

h3 = subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim2), 'FontSize', fontSize);
h3.FontSize = fontSize;

h4 = subplot(3,12,13:17);
scatter(linspace(0,1,M), RxAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_{y|x}', 'FontSize', fontSize);
h4.FontSize = fontSize;

h5 = subplot(3,12,20:24);
scatter(linspace(0,1,M), RyAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_{z|x}', 'FontSize', fontSize);
h5.FontSize = fontSize;

h6 = subplot(3,12,25:29);
scatter(RxAligned,RyAligned); grid on; 
xlabel('r_{y|x}', 'FontSize', fontSize); ylabel('r_{z|x}', 'FontSize', fontSize);  
% title(sprintf('%0.02f/%0.02f/%0.02f/%0.02f/%0.02f', ...
%     ccimVal, pdcorVal, partialCorrVal, cassorVal, hdVal), 'FontSize', fontSize);
title(sprintf('CCIM=%0.02f', ccimVal), 'FontSize', fontSize);
h6.FontSize = fontSize;

h7 = subplot(3,12,32:36);
scatter(pobs(RxAligned),pobs(RyAligned), 'r'); grid on; 
xlabel('F_{r_y}', 'FontSize', fontSize); ylabel('F_{r_z}', 'FontSize', fontSize);
h7.FontSize = fontSize;

%% Do a conditionally dependent test

% Tests Conditional Independence w/ cim and the residuals processing
% algorithm.

clear;
clc;

% rng(6);

M = 500;
noise = 0;
alpha = 0.05;

% Generate data from     Y-->X<--Z
y = rand(M,1);
z = rand(M,1);
% x = y + z + noise*randn(M,1);
% x = (y-0.5).^2 + (z-0.5).^2 + noise*randn(M,1);
% x = sin(4*pi*y)+cos(4*pi*z) + noise*randn(M,1);
% x = nthroot(y,4)+nthroot(z,4) + noise*randn(M,1);
% x = y + (z-0.5).^2 + noise*randn(M,1);
x = (y-0.5).^2 + cos(4*pi*z) + noise*randn(M,1);

cim1 = cim(x, y);
cim2 = cim(x, z);
cim3 = cim(y, z);
% In this graphical model, Y and Z are independent of each other, but when
% conditioned upon X, they become dependent.  Refer to the rules of
% d-separation to see why, this is a V-Structure!
[ccimVal, RxAligned, RyAligned] = ccim(y,z,x);
[~, pdcorVal] = pdcov(y, z, x); pdcorVal = abs(pdcorVal);
partialCorrVal = abs(partialcorr(y,z,x));
data.X = y; data.Y = z; data.Z = x; cassorVal = abs(cassor(data));
hdVal = abs(hd(y, z, x));

fontSize = 20;

h1 = subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', fontSize); ylabel('y', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim1), 'FontSize', fontSize);
h1.FontSize = fontSize;

h2 = subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim3), 'FontSize', fontSize);
h2.FontSize = fontSize;

h3 = subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim2), 'FontSize', fontSize);
h3.FontSize = fontSize;

h4 = subplot(3,12,13:17);
scatter(1:M, RxAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_y', 'FontSize', fontSize);
h4.FontSize = fontSize;

h5 = subplot(3,12,20:24);
scatter(1:M, RyAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_z', 'FontSize', fontSize);
h5.FontSize = fontSize;

h6 = subplot(3,12,25:29);
scatter(RxAligned,RyAligned); grid on; 
xlabel('r_y', 'FontSize', fontSize); ylabel('r_z', 'FontSize', fontSize);  
% title(sprintf('%0.02f/%0.02f/%0.02f/%0.02f/%0.02f', ...
%     ccimVal, pdcor_val, partialCorrVal, cassorVal, hdVal), 'FontSize', fontSize);
title(sprintf('CCIM=%0.02f', ccimVal), 'FontSize', fontSize);
h6.FontSize = fontSize;

subplot(3,12,32:36);
scatter(pobs(RxAligned),pobs(RyAligned), 'r'); grid on; 
xlabel('F_{r_y}', 'FontSize', fontSize); ylabel('F_{r_z}', 'FontSize', fontSize);
h7.FontSize = fontSize;

%% Characterize null distribution {Y indep Z} | X

%% Huawei Plot
clear;
clc;

rng(123);

M = 500;
noise = 0.01;
alpha = 0.05;

% Generate data from     Y<--X-->Z
x = rand(M,1);
% y = x + noise*randn(M,1); z = x + noise*randn(M,1);
% y = 4*(x-0.5).^2 + noise*randn(M,1); z = 4*(x-0.5).^2 + noise*randn(M,1);
% y = sin(4*pi*x) + noise*randn(M,1); z = cos(4*pi*x) + noise*randn(M,1);
% y = nthroot(x, 4) + noise*randn(M,1); z = nthroot(x, 4) + noise*randn(M,1);
% y = x + noise*randn(M,1); z = 4*(x-0.5).^2 + noise*randn(M,1);
y = 4*(x-0.5).^2 + noise*randn(M,1); z = cos(4*pi*x) + noise*randn(M,1);

cim1 = cim(x, y);
cim2 = cim(x, z);
cim3 = cim(y, z);
% ccim conditions X on Y and Z, and sees how related Y and Z
% are to each other ... in this graphical model, they should be UNRELATED
% after the effect of X is removed... i.e. close to independent.  To see
% why, look at the graphical model, Y indep Z | X according to
% d-separation.  So if we condition upon X (i.e. remove teh effect of X on
% Y and Z separately), then we should get independence.
[ccimVal, RxAligned, RyAligned] = ccim(y,z,x);
fontSize = 20;

h1 = subplot(2,2,1);
scatter(x,y); grid on; xlabel('x', 'FontSize', fontSize); ylabel('y', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim1), 'FontSize', fontSize);
h1.FontSize = fontSize;

h2 = subplot(2,2,3);
scatter(y,z); grid on; xlabel('y', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim3), 'FontSize', fontSize);
h2.FontSize = fontSize;

h3 = subplot(2,2,2);
scatter(x,z); grid on; xlabel('x', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('cim=%0.2f', cim2), 'FontSize', fontSize);
h3.FontSize = fontSize;

h4 = subplot(2,2,4);
scatter(RxAligned,RyAligned); grid on; 
xlabel('r_{y|x}', 'FontSize', fontSize); ylabel('r_{z|x}', 'FontSize', fontSize);  
% title(sprintf('%0.02f/%0.02f/%0.02f/%0.02f/%0.02f', ...
%     ccimVal, pdcorVal, partialCorrVal, cassorVal, hdVal), 'FontSize', fontSize);
title(sprintf('CCIM=%0.02f', ccimVal), 'FontSize', fontSize);
h4.FontSize = fontSize;