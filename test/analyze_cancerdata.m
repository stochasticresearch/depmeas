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

% suppress rank-deficient warnings for regression fitting (we aren't
% concerned with that for our current analysis)
spmd
    warning_id =  'MATLAB:rankDeficientMatrix';
    warning('off',warning_id);
end

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
    [~,name,~] = fileparts(fname);
    if(~strcmp(fname,'.') && ~strcmp(fname,'..'))
        fnameWithPath = fullfile(rootDir, 'csv_files', fname);
        fprintf('Processing file=%s\n', fnameWithPath);
        % process this file
        data = csvread(fnameWithPath);
        % transpose the data, b/c computational biologists like to store
        % data in row vectors instead of column vectors ... :/
        data = data';
        [R, RectanglesCell] = pairrsdm( data );
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
        outFile = fullfile(rootDir,'results',[name '.mat']);
        save(outFile, 'R', 'RectanglesCell', 'monotonicityMat', 'uniques', 'numUniques');
    end
end

% restore normal Matlab warning settings
warning('on',warning_id);

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

files = dir(fullfile(rootDir,'results'));
for ii=1:length(files)
    fname = files(ii).name;
    [~,name,~] = fileparts(fname);
    if(~strcmp(fname,'.') && ~strcmp(fname,'..'))
        fnameWithPath = fullfile(rootDir, 'results', fname);
        fprintf('Processing file=%s\n', fnameWithPath);
        load(fnameWithPath);
        
        % figure out what to do with the results
        
        % somehow, need to show the % of monotonicity vs. nonmonotonicity
        % for each dataset, and for the cancer "genre" as a whole also!
    end
end