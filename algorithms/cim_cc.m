function [metric,regionRectangle] = cim_cc(x, y, minScanIncr, alpha)
%CIM - Copula Index for Detecting Dependence and Monotonicity between
%Stochastic Signals.  See associated paper... to be published and preprint
%located here: https://arxiv.org/abs/1703.06686
% Inputs:
%  x - the x variable
%  y - the y variable
%  minscanincr - the minimum scanning increment.  Large
%                values will filter out high frequency dependencies, 
% Outputs:
%  metric - the calculated dependency metric between x and y
%  regionRectagle - the regions detected in the unit square, which each
%                   correspond to a region of monotonicity
% 
%**************************************************************************
%*                                                                        *
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>                *
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

% convert X and Y to pseudo-observations, and scale to be between 0-1
[u,v] = pobs_sorted_cc(x,y);
[v_reverse_sorted,u_reverse_sorted] = pobs_sorted_cc(y,x);

MAX_NUM_RECT = ceil(length(x)/2);
normInvVal = norminv(1-alpha/2);

axisCfgs = [1 2];
ax2minmaxCfgs = { {[0,1]}, {[0,0.5],[0.5,1]} };

% perform a scan pattern while varying U with V full-range, then swap the U-V axes
vecLen = length(axisCfgs)*length(ax2minmaxCfgs);
numScans = ceil(log2(1/minScanIncr))+1;

metricCell = zeros(numScans,MAX_NUM_RECT); 
numPtsCell = zeros(numScans,MAX_NUM_RECT);
rectanglesCell = zeros(numScans,4,MAX_NUM_RECT);
numRectanglesCreatedVec = zeros(1,numScans);

% pre-allocate to max-length for matlab coder speed purposes
metricVecAggr = cell(1,2);
numPtsVecAggr = cell(1,2);
numRectanglesVecAggr = zeros(2,numScans);
rectangleCellAggr = zeros([2 size(rectanglesCell)]);

numMainLoopIter = length(axisCfgs)*length(ax2minmaxCfgs);
maxIICell = -999*ones(numMainLoopIter,2);
rectangleAggr = zeros([numMainLoopIter size(rectangleCellAggr)]);
numRectanglesCreatedMat = -999*ones(numMainLoopIter,2);

% assign cell elements empty stuff to make matlab-coder happy :x
for ii=1:2
    metricVecAggr{ii} = metricCell;
    numPtsVecAggr{ii} = numPtsCell;
end

metrics = zeros(1,vecLen); 
rectangleAggrIdx = 1;
for axisCfg=axisCfgs
    for ax2minmaxCfgsIdx=1:length(ax2minmaxCfgs)
        ax2minmaxCfg = ax2minmaxCfgs{ax2minmaxCfgsIdx};
        
        ax2minmaxCfgLen = length(ax2minmaxCfg);
        
        for ax2mmCfgIdx=1:ax2minmaxCfgLen
            ax2mmCfg = ax2minmaxCfg{ax2mmCfgIdx};
            ax2min = ax2mmCfg(1);
            ax2max = ax2mmCfg(2);

            scanincr = 1;
            for zz=1:numScans
                switch(axisCfg)
                    case 1
                        ax1pts = u; ax2pts = v;
                    otherwise  % changed from case 2 to otherwise for matlab coder
                        ax1pts = v_reverse_sorted; ax2pts = u_reverse_sorted;
                end

                [metricVecTmp, numPtsVecTmp, rectangles, numRectanglesCreated] = ...
                    scanForDep(normInvVal,ax1pts,ax2pts,ax2min,ax2max,scanincr,MAX_NUM_RECT);
                
                metricCell(zz,:) = metricVecTmp;
                numPtsCell(zz,:) = numPtsVecTmp;
                numRectanglesCreatedVec(zz) = numRectanglesCreated;
                rectanglesCell(zz,:,1:size(rectangles,2)) = rectangles;
                
                scanincr = scanincr/2;
            end
            metricVecAggr{ax2mmCfgIdx} = metricCell;
            numPtsVecAggr{ax2mmCfgIdx} = numPtsCell;
            numRectanglesVecAggr(ax2mmCfgIdx,:) = numRectanglesCreatedVec;
            rectangleCellAggr(ax2mmCfgIdx,:,:,:) = rectanglesCell;
        end
        % compute the metric for this.  putting stuff outside the 2nd
        % for-loop allows us to combine the results for {[0,1]} and
        % {[0,0.5][0.5,1]} easily.  metricVecAggr should have the
        % results for when ax2 is {[0,1]} and compute a metric, then it
        % should have the results for {[0,0.5],[0.5,1]} and compute a
        % metric for it.  at the end of processing, we compute a
        % maximum.
        [m,maxIIVec] = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr, numRectanglesVecAggr, ax2minmaxCfgLen);
        maxIICell(rectangleAggrIdx,1:length(maxIIVec)) = maxIIVec;
        metrics(rectangleAggrIdx) = m;
        for ii=1:length(maxIIVec)
            numRectanglesCreatedMat(rectangleAggrIdx,ii) = numRectanglesVecAggr(ii,maxIIVec(ii));
        end
        rectangleAggr(rectangleAggrIdx,:,:,:,:) = rectangleCellAggr;
        
        rectangleAggrIdx = rectangleAggrIdx + 1;
    end
end

[metric,metricMaxIdx] = max(metrics);
if(nargout>1)
    maxII = maxIICell(metricMaxIdx,:);
    if(maxII(2)==-999)
        % means full-scale V configuration was best
        numRects = numRectanglesCreatedMat(metricMaxIdx,1);
        regionRectangle = squeeze(rectangleAggr(metricMaxIdx,1,maxII(1),:,1:numRects));
    else
        % means it was better to break up the rectangle into two halves to
        % maximize dependency.  it is upto the user to merge the rectangles
        % if they deem it necessary
        numRects1 = numRectanglesCreatedMat(metricMaxIdx,1);
        numRects2 = numRectanglesCreatedMat(metricMaxIdx,2);
        regionRectangle1 = squeeze(rectangleAggr(metricMaxIdx,1,maxII(1),:,1:numRects1));
        regionRectangle2 = squeeze(rectangleAggr(metricMaxIdx,2,maxII(2),:,1:numRects2));
        regionRectangle = [regionRectangle1 regionRectangle2];
    end
end

end

function [metric,maxIIVec] = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr, numRectanglesAggr, ax2minmaxCfgLen)
coder.inline('always');
metrics = zeros(2,ax2minmaxCfgLen);
maxIIVec = zeros(1,ax2minmaxCfgLen);

for jj=1:ax2minmaxCfgLen
    groupMetrics = metricVecAggr{jj};
    groupVecLens = numPtsVecAggr{jj};
    groupRectangles = numRectanglesAggr(jj,:);
    
    weightedMetric = -999;
    numPts = -999;
    for ii=1:length(groupRectangles)
        % compute the metric for each group
        numRects = groupRectangles(ii);
        gMetric = groupMetrics(ii,1:numRects);
        gVecLen = groupVecLens(ii,1:numRects);
        
        % compute the weighted metric
        weightedMetricCompute = sum( gMetric.*gVecLen/(sum(gVecLen)) );
        numPtsCompute = sum(gVecLen);
        
        % take the max
        if(weightedMetricCompute>weightedMetric)
            weightedMetric = weightedMetricCompute;
            numPts = numPtsCompute;
            maxIIVec(jj) = ii;
        end
    end
    metrics(1,jj) = weightedMetric;
    metrics(2,jj) = numPts;
end

% combine group metrics
metric = sum( metrics(2,:)/sum(metrics(2,:)).*metrics(1,:) );

end

function [metricVec, numPtsVec, rectangles, rectanglesIdx] = scanForDep(normInvVal, ax1pts, ax2pts, ax2min, ax2max, scanincr, maxNumRect)
%scanForDep - scans for dependencies across the first axis (if you would
%like to scan across the second axis, simply swap the input arguments to 
%this function).
coder.inline('always');

ax1min = 0; ax1max = scanincr;
newRectangle = 1;

metricVec = zeros(1,maxNumRect);
numPtsVec = zeros(1,maxNumRect);
rectangles = zeros(4,maxNumRect);
rectanglesIdx = 1;

metricRectanglePrev = -999;
numPtsPrev = 1;  % should get overwritten
while ax1max<=1
    % find all the points which are contained within this cover rectangle
    matchPts = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max);
    
    numPts = size(matchPts,1);
    if(numPts>=2)   % make sure we have enough points to compute the metric
        % compute the concordance
        tmp_val = taukl_cc( matchPts(:,1),matchPts(:,2),1,0,0);
        taukl_cc_val = 0;  % see: https://www.mathworks.com/help/simulink/ug/calling-matlab-functions.html#bq1h2z9-47
        taukl_cc_val = tmp_val;
        metricRectangle = min(abs(taukl_cc_val),1);  % because we offload to C, sometimes we get
                                                     % values within EPS of 1 
        stdTau = sqrt(4*(1-metricRectangle^2))/sqrt(numPts) * normInvVal;
        if(newRectangle)
            newRectangle = 0;
        else
            % compare to the previous concordance, if there is a change by the
            % threshold amount, rewind the axes of the cover rectangle and 
            if( (metricRectangle < (metricRectanglePrev-stdTau)) )
                metricVec(rectanglesIdx) = metricRectanglePrev;
                numPtsVec(rectanglesIdx) = numPtsPrev;
                rectangles(:,rectanglesIdx) = [ax1min ax1max-scanincr ax2min ax2max];
                rectanglesIdx = rectanglesIdx + 1;
                % start the new cover rectangle
                ax1min = ax1max - scanincr;
                ax1max = ax1min;        % it will be incremented below
                newRectangle = 1;
            end
        end
        metricRectanglePrev = metricRectangle;
        numPtsPrev = numPts;
    end
    ax1max = ax1max + scanincr;
    
    if(ax1max>1)
        if(metricRectanglePrev>=0)
            metricVec(rectanglesIdx) = metricRectanglePrev;
            numPtsVec(rectanglesIdx) = length(ax1pts)-sum(numPtsVec(1:rectanglesIdx));
            rectangles(:,rectanglesIdx) = [ax1min 1 ax2min ax2max];
        end
    end
end

% means we never matched with any points, so compute tau for the range
if(metricRectanglePrev<0)
    tmp_val = taukl_cc( ax1pts,ax2pts,1,0,0 );
    taukl_cc_val = 0;  % see: https://www.mathworks.com/help/simulink/ug/calling-matlab-functions.html#bq1h2z9-47
    taukl_cc_val = tmp_val;
    metricVec(rectanglesIdx) = min(abs(taukl_cc_val),1);
    numPtsVec(rectanglesIdx) = length(ax1pts)-sum(numPtsVec(1:rectanglesIdx));
    rectangles(:,rectanglesIdx) = [0 1 ax2min ax2max];
end

end

function [matchPts] = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max)
coder.inline('always');
matchIdxs = find(ax1pts>ax1min & ax1pts<=ax1max & ax2pts>ax2min & ax2pts<=ax2max);
matchPts = [ax1pts(matchIdxs) ax2pts(matchIdxs)];
end