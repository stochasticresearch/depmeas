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

%% read in the filtered csv file and generate pairwise dep-matrix

clear;
clc;
dbstop if error;

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\climate';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/climate';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/climate';
end
fname = fullfile(rootDir,'normalized_files', 'pollution.csv');

delimiterIn = ',';
headerlinesIn = 1;
pollutionData = importdata(fname,delimiterIn,headerlinesIn);

minSamps = 120;  % ~4 months of data

siteNums = pollutionData.data(:,1);
siteNumsUnique = unique(siteNums);
numSites = length(siteNumsUnique);
pairwiseTimeAlignedData = cell(1,numSites*numSites);
mainIdx = 1;

for ii=1:numSites
    siteI      = siteNumsUnique(ii);
    siteI_idx  = find(pollutionData.data(:,1)==siteI);
    siteIDates = pollutionData.data(siteI_idx,2);
    siteI_NO2  = pollutionData.data(siteI_idx,3);
    siteI_O3   = pollutionData.data(siteI_idx,4);
    siteI_SO2  = pollutionData.data(siteI_idx,5);
    siteI_CO   = pollutionData.data(siteI_idx,6);
    
    iDataChunkIdxs = find(diff(siteIDates)~=1);
    
    for jj=ii+1:numSites
        siteJ = siteNumsUnique(jj);
        siteJ_idx = find(pollutionData.data(:,1)==siteJ);
        siteJDates = pollutionData.data(siteJ_idx,2);
        siteJ_NO2  = pollutionData.data(siteJ_idx,3);
        siteJ_O3   = pollutionData.data(siteJ_idx,4);
        siteJ_SO2  = pollutionData.data(siteJ_idx,5);
        siteJ_CO   = pollutionData.data(siteJ_idx,6);
        
        prevIdx = 1;
        for kk=1:length(iDataChunkIdxs)
            numSampsInChunk = iDataChunkIdxs(kk)-prevIdx+1;
            if(numSampsInChunk>minSamps)
                validIDateCodes = siteIDates(prevIdx+1:iDataChunkIdxs(kk));
                % now try to find corresponding dates that are contiguous
                % in j
                matchVec = ismember(validIDateCodes, siteJDates); matchVec=matchVec';
                % find series of ones that > min # of samples ... this is a
                % valid chunk to process
                is=find(diff([0 matchVec])==1);
                ie=find(diff([matchVec 0])==-1);
                lgt = ie-is+1;
                for ll=1:length(lgt)
                    if(lgt(ll)>minSamps)
                        idxLoI = is(ll);
                        idxHiI = ie(ll);
                        idxLoJ = find(siteJDates==validIDateCodes(idxLoI));
                        idxHiJ = find(siteJDates==validIDateCodes(idxHiI));
                        
                        % perform a cross-check of the datecode alignment
                        % and process the data if it passes the cross-check
                        dateI = validIDateCodes(idxLoI:idxHiI);
                        dateJ = siteJDates(idxLoJ:idxHiJ);
                        crossCheck = sum(dateI-dateJ);

                        if(crossCheck==0)
                            x.datecode = dateI;
                            
                            x.no2_i = siteI_NO2(idxLoI:idxHiI);
                            x.o3_i = siteI_O3(idxLoI:idxHiI);
                            x.so2_i = siteI_SO2(idxLoI:idxHiI);
                            x.co_i = siteI_CO(idxLoI:idxHiI);
                            
                            x.no2_j = siteJ_NO2(idxLoJ:idxHiJ);
                            x.o3_j = siteJ_O3(idxLoJ:idxHiJ);
                            x.so2_j = siteJ_SO2(idxLoJ:idxHiJ);
                            x.co_j = siteJ_CO(idxLoJ:idxHiJ);
                            
                            % store the data for processing
                            pairwiseTimeAlignedData{mainIdx} = x;
                            mainIdx = mainIdx + 1;
                        end
                    end
                end
            end
            prevIdx = iDataChunkIdxs(kk)+1;
        end
    end
end

pairwiseTimeAlignedData = pairwiseTimeAlignedData(~cellfun(@isempty, pairwiseTimeAlignedData));
fname = fullfile(rootDir,'normalized_files', 'pollution_pairwise.mat');
save(fname, 'pairwiseTimeAlignedData','minSamps');

%% process the pairwise data that is time aligned

clear;
clc;
dbstop if error;

warning('off','MATLAB:rankDeficientMatrix')

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\climate';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/climate';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/climate';
end
fname = fullfile(rootDir,'normalized_files', 'pollution_pairwise.mat');
minSamps = 120;  % ~4 months of data
load(fname);

alpha = 0.05;
% https://www.kevinsheppard.com/images/9/95/MFE_Toolbox_Documentation.pdf
% See Page 57 for why the two options below make sense
adfTestType = 0;    % 0 : No deterministic terms
                    % 1 : Constant
                    % 2 : Time Trend
                    % 3 : Constant, DGP assumed to have a time trend
lags = 0;
pairwiseAnalysis = cell(1,length(pairwiseTimeAlignedData));

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the analysis...'),'keepthis','timestamp');
for ii=1:length(pairwiseTimeAlignedData)
    dispstat(sprintf('%0.02f %% complete!',(ii/length(pairwiseTimeAlignedData))*100),'timestamp');
    res.R = zeros(1,28);        % 28 = nchoosek(4,2).  We have 8 different data columns
                            % for each data record, of which we can do a pairwise
                            % analysis on each at a time
    res.RectanglesCell = cell(1,28);
    res.tauklVec = zeros(1,28);
    res.validVec = zeros(1,28);
    res.startDates = zeros(1,28);
    res.stopDates = zeros(1,28);
    res.xData = cell(1,28);
    res.yData = cell(1,28);
    for jj=1:28
        switch(jj)
            case 1
                dataX = pairwiseTimeAlignedData{ii}.no2_i;
                dataY = pairwiseTimeAlignedData{ii}.o3_i;
            case 2
                dataX = pairwiseTimeAlignedData{ii}.no2_i;
                dataY = pairwiseTimeAlignedData{ii}.so2_i;
            case 3
                dataX = pairwiseTimeAlignedData{ii}.no2_i;
                dataY = pairwiseTimeAlignedData{ii}.co_i;
            case 4
                dataX = pairwiseTimeAlignedData{ii}.no2_i;
                dataY = pairwiseTimeAlignedData{ii}.no2_j;
            case 5
                dataX = pairwiseTimeAlignedData{ii}.no2_i;
                dataY = pairwiseTimeAlignedData{ii}.o3_j;
            case 6
                dataX = pairwiseTimeAlignedData{ii}.no2_i;
                dataY = pairwiseTimeAlignedData{ii}.so2_j;
            case 7
                dataX = pairwiseTimeAlignedData{ii}.no2_i;
                dataY = pairwiseTimeAlignedData{ii}.co_j;
            case 8
                dataX = pairwiseTimeAlignedData{ii}.o3_i;
                dataY = pairwiseTimeAlignedData{ii}.so2_i;
            case 9
                dataX = pairwiseTimeAlignedData{ii}.o3_i;
                dataY = pairwiseTimeAlignedData{ii}.co_i;
            case 10
                dataX = pairwiseTimeAlignedData{ii}.o3_i;
                dataY = pairwiseTimeAlignedData{ii}.no2_j;
            case 11
                dataX = pairwiseTimeAlignedData{ii}.o3_i;
                dataY = pairwiseTimeAlignedData{ii}.o3_j;
            case 12
                dataX = pairwiseTimeAlignedData{ii}.o3_i;
                dataY = pairwiseTimeAlignedData{ii}.so2_j;
            case 13
                dataX = pairwiseTimeAlignedData{ii}.o3_i;
                dataY = pairwiseTimeAlignedData{ii}.co_j;
            case 14
                dataX = pairwiseTimeAlignedData{ii}.so2_i;
                dataY = pairwiseTimeAlignedData{ii}.co_i;
            case 15
                dataX = pairwiseTimeAlignedData{ii}.so2_i;
                dataY = pairwiseTimeAlignedData{ii}.no2_j;
            case 16
                dataX = pairwiseTimeAlignedData{ii}.so2_i;
                dataY = pairwiseTimeAlignedData{ii}.o3_j;
            case 17
                dataX = pairwiseTimeAlignedData{ii}.so2_i;
                dataY = pairwiseTimeAlignedData{ii}.so2_j;
            case 18
                dataX = pairwiseTimeAlignedData{ii}.so2_i;
                dataY = pairwiseTimeAlignedData{ii}.co_j;
            case 19
                dataX = pairwiseTimeAlignedData{ii}.co_i;
                dataY = pairwiseTimeAlignedData{ii}.no2_j;
            case 20
                dataX = pairwiseTimeAlignedData{ii}.co_i;
                dataY = pairwiseTimeAlignedData{ii}.o3_j;
            case 21
                dataX = pairwiseTimeAlignedData{ii}.co_i;
                dataY = pairwiseTimeAlignedData{ii}.so2_j;
            case 22
                dataX = pairwiseTimeAlignedData{ii}.co_i;
                dataY = pairwiseTimeAlignedData{ii}.co_j;
            case 23
                dataX = pairwiseTimeAlignedData{ii}.no2_j;
                dataY = pairwiseTimeAlignedData{ii}.o3_j;
            case 24
                dataX = pairwiseTimeAlignedData{ii}.no2_j;
                dataY = pairwiseTimeAlignedData{ii}.so2_j;
            case 25
                dataX = pairwiseTimeAlignedData{ii}.no2_j;
                dataY = pairwiseTimeAlignedData{ii}.co_j;
            case 26
                dataX = pairwiseTimeAlignedData{ii}.o3_j;
                dataY = pairwiseTimeAlignedData{ii}.so2_j;
            case 27
                dataX = pairwiseTimeAlignedData{ii}.o3_j;
                dataY = pairwiseTimeAlignedData{ii}.co_j;
            case 28
                dataX = pairwiseTimeAlignedData{ii}.so2_j;
                dataY = pairwiseTimeAlignedData{ii}.co_j;
        end
        
        % we like column vectors
        % b/c we dropped the invalid data already, we shouldn't have any in
        % the vectors dataX and dataY below...
        dataX = dataX(:);
        dataY = dataY(:);
        if(isempty(find(dataX==-999, 1)) && isempty(find(dataY==-999, 1)))
            % check stationarity of the data
            [~,pvalX] = augdf(dataX,adfTestType,lags);
            [~,pvalY] = augdf(dataY,adfTestType,lags);
            if(pvalX<=alpha && pvalY<=alpha)
                [metric, rectangleCellOut] = cim(dataX,dataY);
                tauklval = taukl(dataX,dataY);
                pval = cimpval(metric, length(dataX));
                if(pval<=alpha)
                    res.R(jj) = metric;
                    res.RectanglesCell{jj} = rectangleCellOut;
                    res.tauklVec(jj) = tauklval;
                    res.startDates(jj) = pairwiseTimeAlignedData{ii}.datecode(1);
                    res.stopDates(jj) = pairwiseTimeAlignedData{ii}.datecode(end);
                    res.validVec(jj) = 1;
                    res.xData{jj} = dataX;
                    res.yData{jj} = dataY;
                end
            end
        end
    end
    pairwiseAnalysis{ii} = res;
end

fname = fullfile(rootDir,'results', 'pollution.mat');
save(fname, 'pairwiseAnalysis');

warning('on','MATLAB:rankDeficientMatrix')

%% Plot the monotonicity results
clear;
clc;
dbstop if error;

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\climate';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/climate';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/climate';
end
fname = fullfile(rootDir,'results', 'pollution.mat');
load(fname);

depThreshVec = [0.01 0.05 0.1 0.15 0.2 0.25];
cimValThresh = 0.4;
finalMonotonicityResults = cell(1,length(depThreshVec));

for zz=1:length(depThreshVec)
    depThresh = depThreshVec(zz);
    fprintf('Processing depThresh=%0.02f\n', depThresh);
    monotonicityResults = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    for ii=1:length(pairwiseAnalysis)
        res = pairwiseAnalysis{ii};
        for jj=1:28
            cimVal = res.R(jj);
            if(res.validVec(jj) && cimVal>=cimValThresh)
                % count the monotonicity after ensuring we didn't overfit
                tauklVal = res.tauklVec(jj);
                percentageDiff = abs(cimVal-tauklVal)/tauklVal;
                if(percentageDiff<=depThresh)
                    numMonotonicRegions = 1;
                else
                    numMonotonicRegions = size(res.RectanglesCell{jj},2);
                end
                % store into the map
                if(isKey(monotonicityResults,numMonotonicRegions))
                    monotonicityResults(numMonotonicRegions) = monotonicityResults(numMonotonicRegions) + 1;
                else
                    monotonicityResults(numMonotonicRegions) = 1;
                end
            end
        end
    end
    keys(monotonicityResults)
    values(monotonicityResults)
    finalMonotonicityResults{zz} = monotonicityResults;
end

save(fullfile(rootDir,'results', 'pollution_finalMonotonicityResults.mat'), 'finalMonotonicityResults', 'depThreshVec');

%% % now plot the results :D

clear;
clc;
dbstop if error;

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\climate';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/climate';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/climate';
end
fname = fullfile(rootDir,'results', 'pollution_finalMonotonicityResults.mat');
load(fname);

maxCount = 0;
for ii=1:length(depThreshVec)
    if(finalMonotonicityResults{ii}.Count>maxCount)
        maxCount = finalMonotonicityResults{ii}.Count;
    end
end

barPlotVec = zeros(maxCount, length(depThreshVec));
for ii=1:length(depThreshVec)
    c = finalMonotonicityResults{ii};
    for jj=1:maxCount
        if(isKey(c,jj))
            val = double(c(jj));    % potentially experiment w/ transforms, 
                                    % log is one
                                    % x/(a+x) is another one
                                    % both diminish the effect of how prevalent
                                    % monotonic dependencies are in the data,
                                    % so for now we leave it as is.
        else
            val = 0;
        end
        valToPlot = val;
        barPlotVec(ii,jj) = valToPlot;
    end
    numTotalPairwiseDepsAnalyzed = sum(barPlotVec(ii,:));
    barPlotVec(ii,:) = barPlotVec(ii,:)/sum(barPlotVec(ii,:)) * 100;
end

barX = depThreshVec*100;
h = bar(barX, barPlotVec,'stacked');
hold on
plot(xlim,[95 95], 'r--', 'LineWidth', 4)
grid on;
fontSize = 20;
xlabel('Tolerance %', 'FontSize', fontSize, 'FontWeight', 'bold')
ylabel('% of total dependencies', 'FontSize', fontSize, 'FontWeight', 'bold');
a = 1:maxCount;
legendCell = cellstr(num2str(a(:)));
legendCell = legendCell.';
lh = legend(legendCell, 'location', 'southeast');
v = get(lh,'title');
set(v,'string',{'# Monotonic', 'Regions'});
title({sprintf('Monotonicity of %d pairwise dependencies analyzed', numTotalPairwiseDepsAnalyzed), ...
       'for Air-Quality Indicators'}, ...
    'FontSize', fontSize, 'FontWeight', 'bold');

xt = get(gca, 'XTick');
set(gca, 'XTick', barX);
set(gca, 'YTick', [20 40 60 80 95]);
ylim([0 100])
set(gca, 'FontSize', 28)