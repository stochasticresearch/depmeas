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

%% Analyze the pair-wise RSDM metrics for as much alingned data as possible
% for all the available stock prices & save results

clear;
clc;
dbstop if error;

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\stocks';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/stocks';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/stocks';
end

% list all available data
files = dir(fullfile(rootDir,'normalized_files'));
numStocksToProcess = length(files)-2;   % the -2 is to exclude '.' and '..'
stocksData = cell(numStocksToProcess,4);   
                                        % {1} = stock name
                                        % {2} = raw data
                                        % {3} = first difference of closing
                                        %       price
                                        % {4} = results of augmented
                                        %       dickey-fuller test of
                                        %       stationarity

% https://www.kevinsheppard.com/images/9/95/MFE_Toolbox_Documentation.pdf
% See Page 57 for why the two options below make sense
adfTestType = 0;    % 0 : No deterministic terms
                    % 1 : Constant
                    % 2 : Time Trend
                    % 3 : Constant, DGP assumed to have a time trend
lags = 0;           % the # of lags to test for
                                        
% read it all in
jj = 1;
for ii=1:length(files)
    fname = files(ii).name;
    [~,name,ext] = fileparts(fname);
    if(strcmpi(ext,'.csv'))
        fnameWithPath = fullfile(rootDir, 'normalized_files', fname);
        fprintf('Processing file=%s\n', fnameWithPath);
        % process this file
        fid = fopen(fnameWithPath);
        header = fgetl(fid);
        data = textscan(fid, '%*s %f %f %f %f %*[^\n]', 'delimiter', ',');
        fclose(fid);
        
        stocksData{jj,1} = name;
        stocksData{jj,2} = data;
        closePrice = data{4};
        % compute returns data
        returnsData = closePrice(1:end-1)-closePrice(2:end);
        stocksData{jj,3} = returnsData;
        % compute augmented dickey-fuller test of stationarity on the
        % returns
        [adftest.stat,adftest.pval,adftest.critval,adftest.resid] = ...
            augdf(returnsData,adfTestType,lags);
        stocksData{jj,4} = adftest;
        
        jj = jj + 1;
    end
end

% start the parallel pool if none exist
p = gcp();

% suppress rank-deficient warnings for regression fitting (we aren't
% concerned with that for our current analysis)
pctRunOnAll warning('off','MATLAB:rankDeficientMatrix');

% now perform monotonicity analysis
R = zeros(numStocksToProcess,numStocksToProcess);
RectanglesCell = cell(numStocksToProcess,numStocksToProcess);
for ii=1:numStocksToProcess
    fprintf('Processing Data # %d\n', ii);
    parfor jj=ii+1:numStocksToProcess
        returns_i = stocksData{ii,3};
        returns_j = stocksData{jj,3};
        numSampsToProcess = min(length(returns_i),length(returns_j));
        returns_i = returns_i(1:numSampsToProcess);
        returns_j = returns_j(1:numSampsToProcess);
        [metric, rectangleCellOut] = rsdm(returns_i,returns_j);
        
        R(ii,jj) = metric; 
        RectanglesCell{ii,jj} = rectangleCellOut;
    end
end

for ii=1:numStocksToProcess
    for jj=ii+1:numStocksToProcess
        R(jj,ii) = R(ii,jj); 
        R(ii,ii) = 1;
        RectanglesCell{jj,ii} = RectanglesCell{ii,jj};
    end
end

sz = size(RectanglesCell);
monotonicityMat = zeros(sz);
for jj=1:sz(1)
    for kk=jj+1:sz(2)
        rco = RectanglesCell{jj,kk};
        monotonicityMat(jj,kk) = size(rco,2);
        monotonicityMat(kk,jj) = size(rco,2);
    end
end
% count the proportion of pairwise dependencies that are monotonic,
% and the # of times the nonmonotonicity occurs.  only count the
% upper-triangle
I = triu(monotonicityMat,1)~=0;
monotonicityVec = monotonicityMat(I);
[uniques, numUniques] = count_unique(monotonicityVec);

outFile = fullfile(rootDir,'stocks_results.mat');
save(outFile);

% restore normal Matlab warning settings
pctRunOnAll warning('on','MATLAB:rankDeficientMatrix');

p = gcp;
delete(p)

%% Plot the Monotonicity results of the above
