function [metric, resid, residAssocIdxs] = ...
    cim_s(x, y, varargin)
%CIM_S - Streaming Copula Index for Detecting Dependence and Monotonicity between
%Stochastic Signals.  See associated paper... to be published and preprint
%located here: 
% Inputs:
%  x - the x variable
%  y - the y variable
%  varargin{1} - minscanincr - the minimum scanning increment.  Large
%                              values will filter out high frequency
%                              dependencies, small values decrease the
%                              statistical power of the dependency metric
%  varargin{2} - diffthresh  - the threshold at which a change in
%                              concordance amount is detected.  Larger
%                              values are more robust to noise, but tend to
%                              miss high frequency changes.
%  varargin{3} - alpha       - the value used to determine significance
%                              level of a box's concordance level
% Outputs:
%  metric - the calculated dependency metric between x and y
%  resid  - the residual between the estimated concordance boxes and the
%           observed statistical variables.  Each concordance box's
%           residuals are provided separately
%  residAssocIdxs - the indices of the independent variable associated with
%                   each residual point, this is used by rscdm for residual
%                   alignment.
%
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

% default values
minscanincr = 0.05;
diffthresh = 120;
alpha = 0.1;

% overwrite defaults with user-inputted values
nVarargin = length(varargin);
switch nVarargin
    case 0
    case 1
        minscanincr = varargin{1};
    case 2
        minscanincr = varargin{1};
        diffthresh = varargin{2};
    otherwise
        % the >=3 case
        minscanincr = varargin{1};
        diffthresh = varargin{2};
        alpha = varargin{3};
end

% convert X and Y to pseudo-observations, and scale to be between 0-1
M = length(x);
u = pobs(x)*(M+1)/M;
v = pobs(y)*(M+1)/M;

axisCfgs = [1 2];
ax2minmaxCfgs = { {[0,1]}, {[0,0.5],[0.5,1]} };

% perform a scan pattern while varying U with V full-range, then swap the U-V axes
metrics = [];
rectangleAggr = {};  rectangleAggrIdx = 1;
for axisCfg=axisCfgs
    for ax2minmaxCfgsIdx=1:length(ax2minmaxCfgs)
        ax2minmaxCfg = ax2minmaxCfgs{ax2minmaxCfgsIdx};
        metricVecAggr = cell(1,length(ax2minmaxCfg));
        numPtsVecAggr = cell(1,length(ax2minmaxCfg));
        rectangleCellAggr = cell(1,length(ax2minmaxCfg));
        for ax2mmCfgIdx=1:length(ax2minmaxCfg)
            ax2mmCfg = ax2minmaxCfg{ax2mmCfgIdx};
            ax2min = ax2mmCfg(1);
            ax2max = ax2mmCfg(2);

            % make the taukl_s object, so that we can efficiently
            % (re)process the data
            switch(axisCfg)
                case 1
                    if(ax2min==0 && ax2max==1)
                        ktsObj = taukl_s(u,v);
                    else
                        I = find(v>ax2min & v<ax2max);
                        uu = u(I); vv = v(I);
                        ktsObj = taukl_s(uu,vv);
                    end
                case 2
                    if(ax2min==0 && ax2max==1)
                        ktsObj = taukl_s(v,u);
                    else
                        I = find(u>ax2min & u<ax2max);
                        uu = u(I); vv = v(I);
                        ktsObj = taukl_s(vv,uu);
                    end
            end
            
            metricCell = {}; numPtsCell = {};
            rectanglesCell = {}; rectanglesCellIdx = 1;
            scanincr = 1;
            while(scanincr>=minscanincr)
                ktsObj.resetState();
                [metricVecTmp, numPtsVecTmp, rectangles] = ...
                    scanForDep(ktsObj,ax2min,ax2max,scanincr,diffthresh,alpha);
                
                metricCell{rectanglesCellIdx} = metricVecTmp;
                numPtsCell{rectanglesCellIdx} = numPtsVecTmp;
                rectanglesCell{rectanglesCellIdx} = rectangles; 
                
                rectanglesCellIdx = rectanglesCellIdx + 1;

                scanincr = scanincr/2;
            end
            metricVecAggr{ax2mmCfgIdx} = metricCell;
            numPtsVecAggr{ax2mmCfgIdx} = numPtsCell;
            rectangleCellAggr{ax2mmCfgIdx} = rectanglesCell;
        end
        % compute the metric for this.  putting stuff outside the 2nd
        % for-loop allows us to combine th results for {[0,1]} and
        % {[0,0.5][0.5,1]} easily.  metricVecAggr should have the
        % results for when ax2 is {[0,1]} and compute a metric, then it
        % should have the results for {[0,0.5],[0.5,1]} and compute a
        % metric for it.  at the end of processing, we compute a
        % maximum.
        m = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr);
        metrics = [metrics m];
        rectangleAggr{rectangleAggrIdx} = rectangleCellAggr;
        rectangleAggrIdx = rectangleAggrIdx + 1;
    end
end

metric = max(metrics);

if(nargout>1)
    ax1pts = u; ax2pts = v;
    orientation = determineOrientation(u,v);
    
    % try to auto-align
    [~,I] = max(metrics(1:length(metrics)/2));
    rectangleCell = rectangleAggr{I};
    
    % if orientation disagree's with the maximum metric, we have to swap
    % the u/v box configurations to align with the rectangle configurations
    % that we will test for
    if(orientation)
        ax1pts = v; ax2pts = u; % flip orientation    
        % swap the rectangles
        for ii=1:length(rectangleCell)
            tmp = rectangleCell{ii};
            for jj=1:length(tmp)
                r = tmp{jj};
                rr = [r(3:4,:); r(1:2,:)];
                rectangleCell{ii}{jj} = rr;
            end
        end
    elseif(~orientation && I>(length(metrics)/2))
        % swap the rectangles
        for ii=1:length(rectangleCell)
            tmp = rectangleCell{ii};
            for jj=1:length(tmp)
                r = tmp{jj};
                rr = [r(3:4,:); r(1:2,:)];
                rectangleCell{ii}{jj} = rr;
            end
        end
    end
    
    % now flatten the axes (we do this without condition b/c we are
    % expecting that we are only residual aligning functional dependencies.
    % how to align non-functional dependencies is unclear at the moment so
    % we don't make any gaurantees in that scenario at the moment!.
    tmp = rectangleCell{1};
    % keep the legacy format... but this is inefficient :(
    rectangleCellNew = {};
    rectangleCellNew{1} = cell(1,length(tmp));
    for jj=1:length(tmp)
        tmpVec = [];
        for kk=1:length(rectangleCell)
            tmpVec = [tmpVec rectangleCell{kk}{jj}];
        end
        rectangleCellNew{1}{jj} = tmpVec;
    end
    rectangleCell = rectangleCellNew;
    
    [resid, residAssocIdxs] = ...
        computeResidual(rectangleCell, ax1pts, ax2pts);
end

end

function [orientation] = determineOrientation(ax1pts, ax2pts, checkBoxSize)
% if this function returns 0, then ax1pts is the independent variable. 
% else, ax2pts is the independent variable.

if(nargin<3)
    checkBoxSize = 0.1;
end

range_ax1_ax2 = zeros(1,round(1/checkBoxSize));
range_ax2_ax1 = zeros(1,round(1/checkBoxSize));

% first check u/v
ax1min = 0; ax1max = checkBoxSize; ax2min = 0; ax2max = 1;
ii = 1;
while(ax1max<1)
    ax1_match = find(ax1pts>ax1min & ax1pts<ax1max);
    ax2_match = find(ax2pts>ax2min & ax2pts<ax2max);
    matchIdxs = intersect(ax1_match,ax2_match);
    ax2matchPts = ax2pts(matchIdxs);
    
    % compute the number of unique points in ax1
    range_ax1_ax2(ii) = quantile(ax2matchPts, 0.75)-quantile(ax2matchPts,0.25);
    
    ii = ii + 1;
    
    ax1min = ax1min + checkBoxSize;
    ax1max = ax1max + checkBoxSize;
end

% now check v/u
tmp = ax1pts; ax1pts = ax2pts; ax2pts = tmp;
ax1min = 0; ax1max = checkBoxSize; ax2min = 0; ax2max = 1;
ii = 1;
while(ax1max<1)
    ax1_match = find(ax1pts>ax1min & ax1pts<ax1max);
    ax2_match = find(ax2pts>ax2min & ax2pts<ax2max);
    matchIdxs = intersect(ax1_match,ax2_match);
    ax2matchPts = ax2pts(matchIdxs);
    
    % compute the number of unique points in ax1
    range_ax2_ax1(ii) = quantile(ax2matchPts, 0.75)-quantile(ax2matchPts,0.25);
    
    ii = ii + 1;
    
    ax1min = ax1min + checkBoxSize;
    ax1max = ax1max + checkBoxSize;
end


if(sum(range_ax1_ax2) <= sum(range_ax2_ax1))
    orientation = 0;
else
    orientation = 1;
end

end

function [resid, residAssocIdxs] = ...
    computeResidual(rectangleCellAggr, ax1pts, ax2pts)

% the cell array that is passed in rectangles for each of the scan-incr's
% we compute the residual for each of these box configurations, and take
% the one that has the lowest residual (i.e. best fit).

resid = {}; residAssocIdxs = {};
residIdx = 1;
for kk=1:length(rectangleCellAggr)
    rectangleCell = rectangleCellAggr{kk};
    residualCell = cell(1,length(rectangleCell));
    residAssocIdxsCell = cell(1,length(rectangleCell));
    for ii=1:length(rectangleCell)
        rectangles = rectangleCell{ii};
        residualMat = cell(1,size(rectangles,2)); 
        residAssocIdxsMat = cell(1,size(rectangles,2));
        residualMatIdx = 1;
        for jj=1:size(rectangles,2)
            rectangle = rectangles(:,jj);
            ax1min = rectangle(1); ax1max = rectangle(2);
            ax2min = rectangle(3); ax2max = rectangle(4);

            ax1_match = find(ax1pts>ax1min & ax1pts<=ax1max);
            ax2_match = find(ax2pts>ax2min & ax2pts<=ax2max);
            matchIdxs = intersect(ax1_match,ax2_match);
            matchPts = [ax1pts(matchIdxs) ax2pts(matchIdxs)];

            XX = [ones(length(matchIdxs),1) matchPts(:,1)];
            zz = XX\matchPts(:,2);
            offset = zz(1); regressionSlope = zz(2);

            yPredict = regressionSlope*matchPts(:,1)+offset;
            residualVec = matchPts(:,2) - yPredict;
            
            residualMat{residualMatIdx} = residualVec;
            residAssocIdxsMat{residualMatIdx} = matchIdxs;
            residualMatIdx = residualMatIdx + 1;
        end
        residualCell{ii} = residualMat;
        residAssocIdxsCell{ii} = residAssocIdxsMat;
    end

    % compute the sum squared of all the residuals, and take the minimum
    minErr = Inf; minErrIdx = 1;
    for ii=1:length(residualCell)
        residualMat = residualCell{ii};
        err = 0;
        for jj=1:length(residualMat)
            err = err + sum(residualMat{jj}.^2);
        end
        if(err<minErr)
            minErr = err;
            minErrIdx = ii;
        end
    end

    resid{residIdx} = residualCell{minErrIdx};
    residAssocIdxs{residIdx} = residAssocIdxsCell{minErrIdx};
    residIdx = residIdx + 1;
end

end

function [metric] = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr)

metrics = zeros(2,length(metricVecAggr));
for ii=1:length(metricVecAggr)
    groupMetrics = metricVecAggr{ii};
    groupVecLens = numPtsVecAggr{ii};
    
    weightedMetric = -999;
    numPts = -999;
    for jj=1:length(groupMetrics)
        % compute the metric for each group
        gMetric = groupMetrics{jj};
        gVecLen = groupVecLens{jj};
        
        % compute the weighted metric
        weightedMetricCompute = sum( gMetric.*gVecLen/(sum(gVecLen)) );
        numPtsCompute = sum(gVecLen);
        
        % take the max
        if(weightedMetricCompute>weightedMetric)
            weightedMetric = weightedMetricCompute;
            numPts = numPtsCompute;
        end
    end
    metrics(1,ii) = weightedMetric;
    metrics(2,ii) = numPts;
end

% combine group metrics
metric = sum( metrics(2,:)/sum(metrics(2,:)).*metrics(1,:) );

end

function [metricVec, numPtsVec, rectangles] = scanForDep(ktsObj, ax2min, ax2max, scanincr, diffthresh, alpha)
%scanForDep - scans for dependencies across the first axis (if you would
%like to scan across the second axis, simply swap the input arguments to 
%this function).

numConsumePerIter = ceil(ktsObj.M/(1.0/scanincr))-1;
if(numConsumePerIter>=ktsObj.M)
    numConsumePerIter = ktsObj.M-1;
end

ax1min = 0; ax1max = scanincr;
newRectangle = 1;
metricVec = [];
numPtsVec = [];
rectangles = []; rectanglesIdx = 1;
numPts = 0;
while ax1max<=1    
    numPts = numPts + numConsumePerIter;
    if(numPts>=2)   % make sure we have enough points to compute the metric
        % compute the concordance
        metricRectangle = abs(ktsObj.consume(numConsumePerIter));
        zsc  = metricRectangle./sqrt( (2*(2*numPts+5))./(9*numPts.*(numPts-1)) );   
        zscDev = abs(zsc-norminv(1-alpha));

        if(newRectangle)
            newRectangle = 0;
        else
            % compare to the previous concordance, if there is a change by the
            % threshold amount, rewind the axes of the cover rectangle and 
            percentageChange = (metricRectangle-metricRectanglePrev)/metricRectanglePrev*100;
            diffthreshAdaptive = (1/zscDev)*diffthresh;
            if(percentageChange<(-1*diffthreshAdaptive))
                metricVec = [metricVec metricRectanglePrev];
                numPtsVec = [numPtsVec numPtsPrev];
                rectangles(:,rectanglesIdx) = [ax1min ax1max-scanincr ax2min ax2max]; rectanglesIdx = rectanglesIdx + 1;
                % start the new cover rectangle
                ax1min = ax1max - scanincr;
                ax1max = ax1min;        % it will be incremented below
                newRectangle = 1; numPts = 0;
                ktsObj.rewind();
                ktsObj.clearState();
            end
        end

        metricRectanglePrev = metricRectangle;
        numPtsPrev = numPts;
    end
    ax1max = ax1max + scanincr;
    
    if(ax1max>1)
        if(exist('metricRectanglePrev','var'))
            metricVec = [metricVec metricRectanglePrev];
            numPtsVec = [numPtsVec numPts];
            rectangles(:,rectanglesIdx) = [ax1min 1 ax2min ax2max]; rectanglesIdx = rectanglesIdx + 1;
        end
    end
end


end