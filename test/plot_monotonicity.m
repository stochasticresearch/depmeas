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

%% Plot the monotonicty results of climate/cancer/stocks in 1 pie-chart 
% that we will insert into the figure for the paper.

clear;
clc;
dbstop if error;

depThresh = 0.05;
fontSize = 36;

%>>>>>>>>>>>>>>>>>>>>>>>> PROCESS CANCER DATA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\cancer';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/cancer';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/cancer';
end
load(fullfile(rootDir,'results', 'finalMonotonicityResults.mat'));
resultsIdx = find(depThresh==depThreshVec);
if(isempty(resultsIdx))
    error('Invalid depThresh chosen!');
end

cancerMonotonicityResults = double(cell2mat(finalMonotonicityResults{resultsIdx}.values()));
% condense the 5 categories into 3 for plot uniformity ...
cancerMonotonicityResults = [cancerMonotonicityResults(1) ...
                             cancerMonotonicityResults(2) ...
                             sum(cancerMonotonicityResults(3:5))];
explode = zeros(1,length(cancerMonotonicityResults)); explode(2) = 1;
%subplot(1,3,1);
figure;
p1 = pie(cancerMonotonicityResults, explode);
p1(2).FontSize = fontSize;
p1(4).FontSize = fontSize;
p1(6).String = '';
fprintf('Num cancer Dependencies Analyzed = %d\n', sum(cancerMonotonicityResults));


%>>>>>>>>>>>>>>>>>>>>>>>> PROCESS STOCKS DATA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\stocks';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/stocks';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/stocks';
end
load(fullfile(rootDir,'finalMonotonicityResults.mat'));
resultsIdx = find(depThresh==depThreshVec);
if(isempty(resultsIdx))
    error('Invalid depThresh chosen!');
end

stocksMonotonicityResults = double(cell2mat(finalMonotonicityResults{resultsIdx}.values()));
explode = zeros(1,length(stocksMonotonicityResults)); explode(2) = 1;
%subplot(1,3,2);
figure;
p2 = pie(stocksMonotonicityResults, explode);
p2(2).FontSize = fontSize;
p2(4).FontSize = fontSize;
p2(6).String = '';
fprintf('Num Stocks Dependencies Analyzed = %d\n', sum(stocksMonotonicityResults));

%>>>>>>>>>>>>>>>>>>>>>>>> PROCESS CLIMATE DATA <<<<<<<<<<<<<<<<<<<<<<<<<<<<
if(ispc)
    rootDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\climate\\results';
elseif(ismac)
    rootDir = '/Users/Kiran/ownCloud/PhD/sim_results/climate/results';
else
    rootDir = '/home/kiran/ownCloud/PhD/sim_results/climate/results';
end
% first get el-nino
load(fullfile(rootDir,'elnino_finalMonotonicityResults.mat'));
resultsIdx = find(depThresh==depThreshVec);
if(isempty(resultsIdx))
    error('Invalid depThresh chosen!');
end
elninoMonotonicityResults = double(cell2mat(finalMonotonicityResults{resultsIdx}.values()));
% now get land temperatures
load(fullfile(rootDir,'landTemperatures_finalMonotonicityResults.mat'));
resultsIdx = find(depThresh==depThreshVec);
if(isempty(resultsIdx))
    error('Invalid depThresh chosen!');
end
landtemperaturesMonotonicityResults = double(cell2mat(finalMonotonicityResults{resultsIdx}.values()));
% now get pollution
load(fullfile(rootDir,'pollution_finalMonotonicityResults.mat'));
resultsIdx = find(depThresh==depThreshVec);
if(isempty(resultsIdx))
    error('Invalid depThresh chosen!');
end
pollutionMonotonicityResults = double(cell2mat(finalMonotonicityResults{resultsIdx}.values()));
% aggregate into climate results
climateResults = elninoMonotonicityResults + landtemperaturesMonotonicityResults + pollutionMonotonicityResults;
explode = zeros(1,length(climateResults)); explode(2) = 1;
figure;
p3 = pie(climateResults, explode);
p3(2).FontSize = fontSize;
p3(4).FontSize = fontSize;
p3(6).String = '';
fprintf('Num Climate Dependencies Analyzed = %d\n', sum(climateResults));