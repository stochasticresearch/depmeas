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

%% RSCDM Merge CD script

clear;
clc;

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
elseif(isunix)
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(rootDir, 'rscdm_CD_gamma_1_25.mat');
f2 = fullfile(rootDir, 'rscdm_CD_gamma_26_50.mat');
f3 = fullfile(rootDir, 'rscdm_CD_gamma_51_75.mat');
f4 = fullfile(rootDir, 'rscdm_CD_gamma_76_101.mat');

load(f2);
cmaResults_f2 = cmaResultsMat;
cmaSurrResults_f2 = cmaSurrResultsMat;
pcorrResults_f2 = pcorrResultsMat;
pcorrPvalResults_f2 = pcorrPValMat;
pdcorrResults_f2 = pdcorrResultsMat;
pdcorrPvalResults_f2 = pdcorrPValMat;
rscdmResults_f2 = rscdmResultsMat;
rsdmResults_f2 = rsdmResultsMat;
clearvars cmaResultsMat cmaSurrResultsMat pcorrResultsMat pcorrPValMat
clearvars pdcorrResultsMat pdcorrPValMat rscdmResultsMat rsdmResultsMat
load(f3);
cmaResults_f3 = cmaResultsMat;
cmaSurrResults_f3 = cmaSurrResultsMat;
pcorrResults_f3 = pcorrResultsMat;
pcorrPvalResults_f3 = pcorrPValMat;
pdcorrResults_f3 = pdcorrResultsMat;
pdcorrPvalResults_f3 = pdcorrPValMat;
rscdmResults_f3 = rscdmResultsMat;
rsdmResults_f3 = rsdmResultsMat;
clearvars cmaResultsMat cmaSurrResultsMat pcorrResultsMat pcorrPValMat
clearvars pdcorrResultsMat pdcorrPValMat rscdmResultsMat rsdmResultsMat
load(f4);
cmaResults_f4 = cmaResultsMat;
cmaSurrResults_f4 = cmaSurrResultsMat;
pcorrResults_f4 = pcorrResultsMat;
pcorrPvalResults_f4 = pcorrPValMat;
pdcorrResults_f4 = pdcorrResultsMat;
pdcorrPvalResults_f4 = pdcorrPValMat;
rscdmResults_f4 = rscdmResultsMat;
rsdmResults_f4 = rsdmResultsMat;
clearvars cmaResultsMat cmaSurrResultsMat pcorrResultsMat pcorrPValMat
clearvars pdcorrResultsMat pdcorrPValMat rscdmResultsMat rsdmResultsMat
load(f1);

gammaVec = 0:.01:1; numDepTypes = 6; nsim = 500;
% merge 26-50
for gammaIdx=26:50
    for jj=1:numDepTypes
        for ii=1:nsim
            cmaResultsMat(gammaIdx,jj,ii) = cmaResults_f2(gammaIdx,jj,ii);
            cmaSurrResultsMat(gammaIdx,jj,ii) = cmaSurrResults_f2(gammaIdx,jj,ii);
            pcorrResultsMat(gammaIdx,jj,ii) = pcorrResults_f2(gammaIdx,jj,ii);
            pcorrPValMat(gammaIdx,jj,ii) = pcorrPvalResults_f2(gammaIdx,jj,ii);
            pdcorrResultsMat(gammaIdx,jj,ii) = pdcorrResults_f2(gammaIdx,jj,ii);
            pdcorrPValMat(gammaIdx,jj,ii) = pdcorrPvalResults_f2(gammaIdx,jj,ii);
            rscdmResultsMat(gammaIdx,jj,ii) = rscdmResults_f2(gammaIdx,jj,ii);
            rsdmResultsMat(gammaIdx,jj,ii) = rsdmResults_f2(gammaIdx,jj,ii);
        end
    end
end

% merge 51-75
for gammaIdx=51:75
    for jj=1:numDepTypes
        for ii=1:nsim
            cmaResultsMat(gammaIdx,jj,ii) = cmaResults_f3(gammaIdx,jj,ii);
            cmaSurrResultsMat(gammaIdx,jj,ii) = cmaSurrResults_f3(gammaIdx,jj,ii);
            pcorrResultsMat(gammaIdx,jj,ii) = pcorrResults_f3(gammaIdx,jj,ii);
            pcorrPValMat(gammaIdx,jj,ii) = pcorrPvalResults_f3(gammaIdx,jj,ii);
            pdcorrResultsMat(gammaIdx,jj,ii) = pdcorrResults_f3(gammaIdx,jj,ii);
            pdcorrPValMat(gammaIdx,jj,ii) = pdcorrPvalResults_f3(gammaIdx,jj,ii);
            rscdmResultsMat(gammaIdx,jj,ii) = rscdmResults_f3(gammaIdx,jj,ii);
            rsdmResultsMat(gammaIdx,jj,ii) = rsdmResults_f3(gammaIdx,jj,ii);
        end
    end
end

% merge 76-101
for gammaIdx=76:101
    for jj=1:numDepTypes
        for ii=1:nsim
            cmaResultsMat(gammaIdx,jj,ii) = cmaResults_f4(gammaIdx,jj,ii);
            cmaSurrResultsMat(gammaIdx,jj,ii) = cmaSurrResults_f4(gammaIdx,jj,ii);
            pcorrResultsMat(gammaIdx,jj,ii) = pcorrResults_f4(gammaIdx,jj,ii);
            pcorrPValMat(gammaIdx,jj,ii) = pcorrPvalResults_f4(gammaIdx,jj,ii);
            pdcorrResultsMat(gammaIdx,jj,ii) = pdcorrResults_f4(gammaIdx,jj,ii);
            pdcorrPValMat(gammaIdx,jj,ii) = pdcorrPvalResults_f4(gammaIdx,jj,ii);
            rscdmResultsMat(gammaIdx,jj,ii) = rscdmResults_f4(gammaIdx,jj,ii);
            rsdmResultsMat(gammaIdx,jj,ii) = rsdmResults_f4(gammaIdx,jj,ii);
        end
    end
end

fout = fullfile(rootDir, 'rscdm_CD_new.mat');
significance_thresh = 0.05;
save(fout, 'cmaResultsMat', 'cmaSurrResultsMat', 'pcorrResultsMat', 'pcorrPValMat', ...
    'pdcorrResultsMat', 'pdcorrPValMat', 'rscdmResultsMat', 'rsdmResultsMat', ...
    'gammaVec', 'numDepTypes', 'nsim', 'significance_thresh');