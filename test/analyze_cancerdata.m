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

%%
% Compute monotonicity results of cancer datasets provided by the 
% Broad Institute

clear;
clc;
dbstop if error;

parpool;
% suppress rank-deficient warnings for regression fitting (we aren't
% concerned with that for our current analysis)
pctRunOnAll warning('off','MATLAB:rankDeficientMatrix');

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\cancer';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/cancer';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/cancer';
end

files = dir(fullfile(rootDir,'csv_files'));
for ii=1:length(files)
    fname = files(ii).name;
    [~,name,ext] = fileparts(fname);
    if(strcmpi(ext,'.csv'))
        fnameWithPath = fullfile(rootDir, 'csv_files', fname);
        outFile = fullfile(rootDir,'results',[name '.mat']);
        if(exist(outFile,'file'))
            % we have already processed the results, no need to recompute
            continue;
        end
        fprintf('Processing file=%s\n', fnameWithPath);
        % process this file
        data = csvread(fnameWithPath);
        % transpose the data, b/c computational biologists like to store
        % data in row vectors instead of column vectors ... :/
        data = data';
        if(size(data,2)>100)
            data = data(:,1:100);
        end
        [R, RectanglesCell] = paircim( data );
        % compute the "monotonicity" of the data
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
        save(outFile, 'R', 'RectanglesCell', 'monotonicityMat', 'uniques', 'numUniques', 'data');
    end
end

% restore normal Matlab warning settings
pctRunOnAll warning('on','MATLAB:rankDeficientMatrix');

p = gcp;
delete(p)

%%
% analyze the results of monotonicity

clear;
clc;

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\cancer';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/cancer';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/cancer';
end

alpha = 0.05;       % for p-value comparison
depThreshVec = [0.01 0.05 0.1 0.15 0.2 0.25];
finalMonotonicityResults = cell(1,length(depThreshVec));

for zz=1:length(depThreshVec)
    depThresh = depThreshVec(zz);   % vary the percentage difference tolerated
    
    fprintf('Processing depThresh=%0.02f\n', depThresh);
    
    files = dir(fullfile(rootDir,'results'));
    % first delete any postProcessed files we may already have
    for ii=1:length(files)
        fname = files(ii).name;
        [~,name,ext] = fileparts(fname);
        if(~isempty(strfind(name,'postProcessed')))
            delete(fname);
        end
    end
    % now post-process again with the chosen depThresh
    files = dir(fullfile(rootDir,'results'));
    for ii=1:length(files)
        fname = files(ii).name;
        [~,name,ext] = fileparts(fname);

        % delete any old 'postProcessed files', since we will do the
        % postProcessing again here
        if(strcmpi(ext,'.mat'))
            fnameWithPath = fullfile(rootDir, 'results', fname);
            %fprintf('Processing file=%s\n', fnameWithPath);
            load(fnameWithPath);

            % Find all the pairwise dependencies that were flagged as
            % nonmonotonic, and compare the value of CIM to taukl, if they are
            % sufficiently close, it means that we have overfit, so we can
            % conclude that the dependency is indeed actually monotonic
            nonMonotonicIdx = find(monotonicityMat>1); nonMonotonicIdx = nonMonotonicIdx';
            for idx=nonMonotonicIdx
                [iIdx,jIdx] = ind2sub(size(monotonicityMat), idx);
                % get the CIM value
                cimVal = R(iIdx,jIdx);
                tauklVal = abs(taukl(data(:,iIdx),data(:,jIdx)));
                percentageDiff = abs(cimVal-tauklVal)/tauklVal;
                if(percentageDiff<=depThresh)
                    % means we overfit, and we correct for that here
                    monotonicityMat(iIdx,jIdx) = 1;
                    monotonicityMat(jIdx,iIdx) = 1;
                end
            end
            
            numDataPoints = size(data,1);
            validDepMat = zeros(size(monotonicityMat));
            for jj=1:size(R,1)
                for kk=jj+1:size(R,1)
                    pval = cimpval(R(jj,kk), numDataPoints);
                    if(pval<alpha)
                        validDepMat(jj,kk) = 1; % we flag this as non-independent, and
                                                % thus we will process it
                        validDepMat(kk,jj) = 1;
                    end
                end
            end

            % again, count the # of unique monotonicity levels we have
            I1 = find(triu(monotonicityMat,1)~=0);
            I2 = find(validDepMat~=0);
            intersectI = intersect(I1,I2);
            monotonicityVec = monotonicityMat(intersectI);
            [uniques, numUniques] = count_unique(monotonicityVec);
            numDepsAnalyzed = length(intersectI);

            outFile = fullfile(rootDir, 'results', [name '_postProcessed.mat']);
            save(outFile, 'depThresh', 'R', 'RectanglesCell', 'monotonicityVec', 'uniques', 'numUniques', 'data', 'validDepMat', 'numDepsAnalyzed');
        end
    end

    minNumSamples = 50;
    numFilesAnalyzed = 0;
    % now plot the aggregate monotonicity results
    % now post-process again with the chosen depThresh
    files = dir(fullfile(rootDir,'results'));
    monotonicityResults = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    numTotalDepsAnalyzed = 0;
    for ii=1:length(files)
        fname = files(ii).name;
        [~,name,ext] = fileparts(fname);
        if(~isempty(strfind(name,'postProcessed')) && strcmpi(ext,'.mat'))
            fnameWithPath = fullfile(rootDir, 'results', fname);
            %fprintf('Processing file=%s\n', fnameWithPath);
            load(fnameWithPath);

            % if we have processed > minNumSamples samples at minimum, 
            % aggregate the results
            if(size(data,1)>minNumSamples)
                numFilesAnalyzed = numFilesAnalyzed + 1;
                numTotalDepsAnalyzed = numTotalDepsAnalyzed + numDepsAnalyzed;
                lenUniques = length(uniques);
                for jj=1:lenUniques
                    if(isKey(monotonicityResults, uniques(jj)))
                        monotonicityResults(uniques(jj)) = monotonicityResults(uniques(jj)) + numUniques(jj);
                    else
                        monotonicityResults(uniques(jj)) = numUniques(jj);
                    end
                end
            end
        end
    end
    numTotalDepsAnalyzed
    keys(monotonicityResults)
    values(monotonicityResults)
    finalMonotonicityResults{zz} = monotonicityResults;
end

% save off the results
save(fullfile(rootDir,'results', 'finalMonotonicityResults.mat'), 'finalMonotonicityResults', 'depThreshVec', 'numTotalDepsAnalyzed');

%% plot the final monotonicityResults

clear;
clc;

if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\cancer';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/cancer';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/cancer';
end
load(fullfile(rootDir,'results', 'finalMonotonicityResults.mat'));

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
        if(~isempty(find(cell2mat(c.keys())==jj)))
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
title({sprintf('Monotonicity of %d pairwise dependencies analyzed', ...
    numTotalPairwiseDepsAnalyzed), 'for various Gene-Expression datasets'}, ...
    'FontSize', fontSize, 'FontWeight', 'bold');

xt = get(gca, 'XTick');
set(gca, 'XTick', barX);
set(gca, 'YTick', [20 40 60 80 95]);
ylim([0 100])
set(gca, 'FontSize', 28)