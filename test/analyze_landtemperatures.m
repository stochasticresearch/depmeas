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

%% process the pairwise dependencies for average land temperatures

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
fname = fullfile(rootDir,'normalized_files', 'landTemperatures.mat');

warning('off','MATLAB:rankDeficientMatrix')

countries = who('-file', fname);
numCountries = length(countries);
load(fname);        % load the data

alpha = 0.05;
% https://www.kevinsheppard.com/images/9/95/MFE_Toolbox_Documentation.pdf
% See Page 57 for why the two options below make sense
adfTestType = 0;    % 0 : No deterministic terms
                    % 1 : Constant
                    % 2 : Time Trend
                    % 3 : Constant, DGP assumed to have a time trend
lags = 0;   
% see if data is stationary
validCountriesVec = zeros(1,numCountries);
for ii=1:numCountries
    countryIData = eval(countries{ii});
    countryIData = countryIData(:);
    % find the last occurance of -999.0, which was inserted if no data was
    % found for a certain date by the raw-data processor
    I = find(countryIData==-999.0, 1, 'last' );
    if(~isempty(I))
        countryIData = countryIData(I:end);
    end
    [adftest.stat,adftest.pval,adftest.critval,adftest.resid] = ...
        augdf(countryIData,adfTestType,lags);
    if(adftest.pval<alpha)
        validCountriesVec(ii) = 1;
    end
end

% compute the pairwise CIM metrics for appropriate data
R = zeros(numCountries);
RectanglesCell = cell(numCountries);
tauklMat = zeros(numCountries);
validMat = zeros(numCountries);
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the analysis...'),'keepthis','timestamp');
for ii=1:numCountries
    for jj=ii+1:numCountries  % we can't do parfor here b/c of eval :(
        percentageComplete = jj/numCountries*100;
        dispstat(sprintf('Processing %s --> %0.02f %%',countries{ii},percentageComplete),'timestamp');
        if(validCountriesVec(ii) && validCountriesVec(jj))
            countryIData = eval(countries{ii});
            countryJData = eval(countries{jj});
            countryIData = countryIData(:);
            countryJData = countryJData(:);
            Ii = find(countryIData==-999.0, 1, 'last' );
            if(~isempty(Ii))
                countryIData = countryIData(Ii:end);
            end
            Ij = find(countryJData==-999.0, 1, 'last' );
            if(~isempty(Ij))
                countryJData = countryJData(Ij:end);
            end
            numSampsToProc = min(length(countryIData),length(countryJData));
            countryIData = countryIData(end-numSampsToProc+1:end);
            countryJData = countryJData(end-numSampsToProc+1:end);
            [metric, rectangleCellOut] = cim(countryIData,countryJData);
            pval = cimpval(metric, numSampsToProc);
            if(pval<alpha)
                R(ii,jj) = metric; R(jj,ii) = metric;
                RectanglesCell{ii,jj} = rectangleCellOut; RectanglesCell{jj,ii} = rectangleCellOut;
                validMat(ii,jj) = 1; validMat(jj,ii) = 1;
                tauklVal = taukl(countryIData,countryJData);
                tauklMat(ii,jj) = tauklVal;
            end
        end
    end
    dispstat(sprintf('%s complete!',countries{ii}),'keepthis', 'timestamp');
end

% save the data for post-processing
fname = fullfile(rootDir,'results', 'landTemperaturesResults.mat');
save(fname, 'R', 'RectanglesCell', 'validMat', 'tauklMat');

warning('on','MATLAB:rankDeficientMatrix')
%% Post-Process & plot the monotonicity results

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
fname = fullfile(rootDir,'results', 'landTemperaturesResults.mat');
load(fname);

numCountries = size(R,1);
depThreshVec = [0.01 0.05 0.1 0.15 0.2 0.25];
cimValThresh = 0.4;
finalMonotonicityResults = cell(1,length(depThreshVec));

for zz=1:length(depThreshVec)
    depThresh = depThreshVec(zz);
    fprintf('Processing depThresh=%0.02f\n', depThresh);
    monotonicityResults = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    for ii=1:numCountries
        for jj=ii+1:numCountries
            cimVal = R(ii,jj);
            if(validMat(ii,jj) && cimVal>=cimValThresh)
                % means both data are stationary, and this dependnecy is
                % significant
                
                % ensure we didn't overfit and compute the monotonicity
                tauklVal = tauklMat(ii,jj);
                percentageDiff = abs(cimVal-tauklVal)/tauklVal;
                if(percentageDiff<=depThresh)
                    numMonotonicRegions = 1;
                else
                    numMonotonicRegions = size(RectanglesCell{ii,jj},2);
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

save(fullfile(rootDir,'results', 'landTemperatures_finalMonotonicityResults.mat'), 'numCountries', 'finalMonotonicityResults', 'depThreshVec');

%% 
% now plot the results :D
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
fname = fullfile(rootDir,'results', 'landTemperatures_finalMonotonicityResults.mat');
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
       sprintf('for Average Land Temperatures of %d countries', numCountries)}, ...
    'FontSize', fontSize, 'FontWeight', 'bold');

xt = get(gca, 'XTick');
set(gca, 'XTick', barX);
set(gca, 'YTick', [20 40 60 80 95]);
ylim([0 100])
set(gca, 'FontSize', 28)

%% Manually consider the monotonicity results

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
fname = fullfile(rootDir,'results', 'landTemperaturesResults.mat');
load(fname);
fname = fullfile(rootDir,'normalized_files', 'landTemperatures.mat');
countries = who('-file', fname);
load(fname);

numCountries = size(R,1);
depThreshVec = [0.01 0.05 0.1 0.15 0.2 0.25];
cimValThresh = 0.4;
finalMonotonicityResults = cell(1,length(depThreshVec));

depThresh = depThreshVec(4);
for ii=1:numCountries
    for jj=ii+1:numCountries
        cimVal = R(ii,jj);
        if(validMat(ii,jj) && cimVal>=cimValThresh)
            % means both data are stationary, and this dependnecy is
            % significant

            % ensure we didn't overfit and compute the monotonicity
            tauklVal = tauklMat(ii,jj);
            percentageDiff = abs(cimVal-tauklVal)/tauklVal;
            if(percentageDiff<=depThresh)
                numMonotonicRegions = 1;
            else
                numMonotonicRegions = size(RectanglesCell{ii,jj},2);
                
                countryIData = eval(countries{ii});
                countryJData = eval(countries{jj});
                countryIData = countryIData(:);
                countryJData = countryJData(:);
                Ii = find(countryIData==-999.0, 1, 'last' );
                if(~isempty(Ii))
                    countryIData = countryIData(Ii:end);
                end
                Ij = find(countryJData==-999.0, 1, 'last' );
                if(~isempty(Ij))
                    countryJData = countryJData(Ij:end);
                end
                numSampsToProc = min(length(countryIData),length(countryJData));
                countryIData = countryIData(end-numSampsToProc+1:end);
                countryJData = countryJData(end-numSampsToProc+1:end);
                [metric, rectangleCellOut] = cim(countryIData,countryJData);
                rectangleCellOut
                scatter(pobs(countryIData),pobs(countryJData));
                xlabel(countries{ii}); ylabel(countries{jj});
                pause;
                
            end

        end
    end
end