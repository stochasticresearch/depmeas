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

%% process the pairwise dependencies for el-nino indicators

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
fname = fullfile(rootDir,'normalized_files', 'elnino.mat');

buoys = who('-file',fname);
numBuoys = length(buoys);
load(fname);        % load the data

minSamps = 100;

% matlab automatically grows the cell array size as needed, although this
% is slower
pairwiseTimeAlignedData = cell(1,numBuoys*numBuoys);
mainIdx = 1;
for ii=1:numBuoys
    buoyIData = eval(buoys{ii});
    iDataChunkIdxs = find(diff(buoyIData.datecode)~=1);
    for jj=ii+1:numBuoys
        buoyJData = eval(buoys{jj});
        
        % find overlapping regions of contiguous data that have > minimum
        % samples, 
        prevIdx = 1;
        for kk=1:length(iDataChunkIdxs)
            numSampsInChunk = iDataChunkIdxs(kk)-prevIdx+1;
            if(numSampsInChunk>minSamps)
                validIDateCodes = buoyIData.datecode(prevIdx+1:iDataChunkIdxs(kk));
                % now try to find corresponding dates that are contiguous
                % in j
                matchVec = ismember(validIDateCodes, buoyJData.datecode);
                % find series of ones that > min # of samples ... this is a
                % valid chunk to process
                is=find(diff([0 matchVec])==1);
                ie=find(diff([matchVec 0])==-1);
                lgt = ie-is+1;
                for ll=1:length(lgt)
                    if(lgt(ll)>minSamps)
                        idxLoI = is(ll);
                        idxHiI = ie(ll);
                        idxLoJ = find(buoyJData.datecode==validIDateCodes(idxLoI));
                        idxHiJ = find(buoyJData.datecode==validIDateCodes(idxHiI));
                        
                        % perform a cross-check of the datecode alignment
                        % and process the data if it passes the cross-check
                        dateI = validIDateCodes(idxLoI:idxHiI);
                        dateJ = buoyJData.datecode(idxLoJ:idxHiJ);
                        crossCheck = sum(dateI-dateJ);
%                         fprintf('buoyI[%d,%d]=%d buoyJ[%d,%d]=%d crossCheckSum=%d\n', ...
%                             ii, idxLoI, idxHiI, jj, idxLoJ, idxHiJ, crossCheck);
                        
                        if(crossCheck==0)
                            % means all datecodes align so we store the
                            % data
                            x.datecode = dateI;
                            
                            x.zonalWindsI = buoyIData.zonalWinds(idxLoI:idxHiI);
                            x.airTempI = buoyIData.airTemp(idxLoI:idxHiI);
                            x.humidityI = buoyIData.airTemp(idxLoI:idxHiI);
                            x.meridionalWindsI = buoyIData.meridionalWinds(idxLoI:idxHiI);
                            x.seaSurfaceTempI = buoyIData.seaSurfaceTemp(idxLoI:idxHiI);
                            
                            x.zonalWindsJ = buoyJData.zonalWinds(idxLoJ:idxHiJ);
                            x.airTempJ = buoyJData.airTemp(idxLoJ:idxHiJ);
                            x.humidityJ = buoyJData.airTemp(idxLoJ:idxHiJ);
                            x.meridionalWindsJ = buoyJData.meridionalWinds(idxLoJ:idxHiJ);
                            x.seaSurfaceTempJ = buoyJData.seaSurfaceTemp(idxLoJ:idxHiJ);
                            
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
fname = fullfile(rootDir,'normalized_files', 'elnino_pairwise.mat');
save(fname, 'pairwiseTimeAlignedData','minSamps');
%% process the pairwise data that is available

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
fname = fullfile(rootDir,'normalized_files', 'elnino_pairwise.mat');
minSamps = 100;
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

warning('off','MATLAB:rankDeficientMatrix')

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the analysis...'),'keepthis','timestamp');
for ii=1:length(pairwiseTimeAlignedData)
    dispstat(sprintf('%0.02f %% complete!',(ii/length(pairwiseTimeAlignedData))*100),'timestamp');
    res.R = zeros(1,45);        % 45 = nchoosek(10,2).  We have 10 different data records
                            % for each buoy, of which we can do a pairwise
                            % analysis on each at a time
    res.RectanglesCell = cell(1,45);
    res.tauklVec = zeros(1,45);
    res.validVec = zeros(1,45);
    res.startDates = zeros(1,45);
    res.stopDates = zeros(1,45);
    for jj=1:45
        switch(jj)
            case 1
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.airTempI;
            case 2
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.humidityI;
            case 3
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsI;
            case 4
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
            case 5
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.airTempJ;
            case 6
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.zonalWindsJ;
            case 7
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.humidityJ;
            case 8
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 9
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 10
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.humidityI;
            case 11
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsI;
            case 12
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
            case 13
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.airTempJ;
            case 14
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.zonalWindsJ;
            case 15
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.humidityJ;
            case 16
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 17
                dataX = pairwiseTimeAlignedData{ii}.airTempI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 18
                dataX = pairwiseTimeAlignedData{ii}.humidityI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsI;
            case 19
                dataX = pairwiseTimeAlignedData{ii}.humidityI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
            case 20
                dataX = pairwiseTimeAlignedData{ii}.humidityI;
                dataY = pairwiseTimeAlignedData{ii}.zonalWindsJ;
            case 21
                dataX = pairwiseTimeAlignedData{ii}.humidityI;
                dataY = pairwiseTimeAlignedData{ii}.airTempJ;
            case 22
                dataX = pairwiseTimeAlignedData{ii}.humidityI;
                dataY = pairwiseTimeAlignedData{ii}.humidityJ;
            case 23
                dataX = pairwiseTimeAlignedData{ii}.humidityI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 24
                dataX = pairwiseTimeAlignedData{ii}.humidityI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 25
                dataX = pairwiseTimeAlignedData{ii}.meridionalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
            case 26
                dataX = pairwiseTimeAlignedData{ii}.meridionalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.zonalWindsJ;
            case 27
                dataX = pairwiseTimeAlignedData{ii}.meridionalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.airTempJ;
            case 28
                dataX = pairwiseTimeAlignedData{ii}.meridionalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.humidityJ;
            case 29
                dataX = pairwiseTimeAlignedData{ii}.meridionalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 30
                dataX = pairwiseTimeAlignedData{ii}.meridionalWindsI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 31
                dataX = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
                dataY = pairwiseTimeAlignedData{ii}.zonalWindsJ;
            case 32
                dataX = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
                dataY = pairwiseTimeAlignedData{ii}.airTempJ;
            case 33
                dataX = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
                dataY = pairwiseTimeAlignedData{ii}.humidityJ;
            case 34
                dataX = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 35
                dataX = pairwiseTimeAlignedData{ii}.seaSurfaceTempI;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 36
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsJ;
                dataY = pairwiseTimeAlignedData{ii}.airTempJ;
            case 37
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsJ;
                dataY = pairwiseTimeAlignedData{ii}.humidityJ;
            case 38
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsJ;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 39
                dataX = pairwiseTimeAlignedData{ii}.zonalWindsJ;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 40
                dataX = pairwiseTimeAlignedData{ii}.airTempJ;
                dataY = pairwiseTimeAlignedData{ii}.humidityJ;
            case 41
                dataX = pairwiseTimeAlignedData{ii}.airTempJ;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 42
                dataX = pairwiseTimeAlignedData{ii}.airTempJ;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 43
                dataX = pairwiseTimeAlignedData{ii}.humidityJ;
                dataY = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
            case 44
                dataX = pairwiseTimeAlignedData{ii}.humidityJ;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
            case 45
                dataX = pairwiseTimeAlignedData{ii}.meridionalWindsJ;
                dataY = pairwiseTimeAlignedData{ii}.seaSurfaceTempJ;
        end
        
        % we like column vectors
        dataX = dataX(:);
        dataY = dataY(:);
        
        % get the largest subset of dataX and dataY that is not 
        % missing data, as noted by -999.0 records
        invalidX = find(dataX==-999.0);
        invalidY = find(dataY==-999.0);
        datesX = pairwiseTimeAlignedData{ii}.datecode;
        datesX(invalidX) = Inf;
        datesY = pairwiseTimeAlignedData{ii}.datecode;
        datesY(invalidY) = -Inf;
        mutuallyGoodDates = intersect(datesX,datesY);
        diffVec = diff(mutuallyGoodDates);
        diffI = find(diffVec~=1);       % this is where we have breaks in the data
        
        % take the largest subset
        maxRun = 1;
        maxStartIdx = 1; maxEndIdx = 0;
        prevIdx = 1;
        for kk=1:length(diffI)
            currRunLen = diffI(kk)-prevIdx;
            if(currRunLen>maxRun)
                maxRun = currRunLen;
                maxStartIdx = prevIdx;
                maxEndIdx = diffI(kk);
            end
            prevIdx = diffI(kk);
        end
        
        % if that subset > minSamps, process and store
        if(maxRun>minSamps)
            dataX = dataX(maxStartIdx:maxEndIdx);
            dataY = dataY(maxStartIdx:maxEndIdx);
            
            % check stationarity of the data!
            [~,pvalX] = augdf(dataX,adfTestType,lags);
            [~,pvalY] = augdf(dataY,adfTestType,lags);

            if(pvalX<alpha && pvalY<alpha)
                [metric, rectangleCellOut] = rsdm(dataX,dataY);
                tauklval = taukl(dataX,dataY);
                pval = rsdmpval(metric, length(dataX));
                if(pval<alpha)
                    res.R(jj) = metric;
                    res.RectanglesCell{jj} = rectangleCellOut;
                    res.tauklVec(jj) = tauklval;
                    res.startDates(jj) = datesX(maxStartIdx);
                    res.stopDates(jj) = datesX(maxEndIdx);
                    res.validVec(jj) = 1;
                end
            end
        end
    end
    pairwiseAnalysis{ii} = res;
end

fname = fullfile(rootDir,'results', 'elnino.mat');
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
fname = fullfile(rootDir,'results', 'elnino.mat');
load(fname);

depThreshVec = [0.01 0.05 0.1 0.15 0.2 0.25];
finalMonotonicityResults = cell(1,length(depThreshVec));


for zz=1:length(depThreshVec)
    depThresh = depThreshVec(zz);
    fprintf('Processing depThresh=%0.02f\n', depThresh);
    monotonicityResults = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    for ii=1:length(pairwiseAnalysis)
        res = pairwiseAnalysis{ii};
        for jj=1:45
            if(res.validVec(jj))
                % count the monotonicity after ensuring we didn't overfit
                rsdmVal = res.R(jj);
                tauklVal = res.tauklVec(jj);
                percentageDiff = abs(rsdmVal-tauklVal)/tauklVal;
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

save(fullfile(rootDir,'results', 'elnino_finalMonotonicityResults.mat'), 'finalMonotonicityResults', 'depThreshVec');

%% Plot
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
fname = fullfile(rootDir,'results', 'elnino_finalMonotonicityResults.mat');
load(fname);

% now plot the results :D
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
       'for El-Nino Indicators'}, ...
    'FontSize', fontSize, 'FontWeight', 'bold');

xt = get(gca, 'XTick');
set(gca, 'XTick', barX);
set(gca, 'YTick', [20 40 60 80 95]);
ylim([0 100])
set(gca, 'FontSize', 28)