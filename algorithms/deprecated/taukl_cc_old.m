function [ tau ] = taukl_cc_old( U, V )
%TAUKL - computes a rescaled version of Kendall's tau that preserves
%         the definition of Kendall's tau, but assures that in the 
%         scenario of perfect concordance or discordance for discrete
%         or hybrid datatypes, taucj achieves +/- 1 respectively
% Inputs:
%  U - first variable input.
%  V - second variable input.
%  THE CC VERSION REQUIRES U & V to be sorted by U!
% Outputs:
%  tau - the rescaled version of Kendall's tau
%  
% NOTE: See Theorem 12 of the paper "On Rank Correlation Measures for
%       Non-Continuous Random Variables" - Neslehova (2007) for more info.
%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%* 
%**************************************************************************
len = length(U);

% compute the numerator the tau_hat
% K = 0;
% for k = 1:len-1
%     K = K + sum( sign(U(k)-U(k+1:len)) .* sign(V(k)-V(k+1:len)) );
% end
K = kendallsTauNumer1(U,V);
% K = kendallsTauNumer2(U,V);

if(K==0)
    tau = 0;
    return;
end

% compute the denominator ... compute the # of unique values of U and V and
% how many times each of those unique values occur
% uniqueUValues = unique(U);
% uniqueVValues = unique(V);
[uniqueUValues, uniqueUCounts] = uniqueSorted(U);
[uniqueVValues, uniqueVCounts] = uniqueSorted(V);

if((length(uniqueUValues) >= len/2) && (length(uniqueVValues) >= len/2))
    % means we can reasonably assume that we have continuous data
    % for data-sizes >=50
    tau = K/(len*(len-1)/2);
    return;
end

% % uniqueUCounts = zeros(1,length(uniqueU));
% % uniqueVCounts = zeros(1,length(uniqueV));
% % 
% % % TODO: we can combine the loops below after verification of functionality
% % for ii=1:length(uniqueU)
% %     uniqueUCounts(ii) = sum(U==uniqueU(ii));
% % end
% % 
% % for ii=1:length(uniqueV)
% %     uniqueVCounts(ii) = sum(V==uniqueV(ii));
% % end

u = 0;
k = 2;
for ii=1:length(uniqueUValues)
    n = uniqueUCounts(ii);
    if(k<=n)
        addVal = nchoosek(n,k);
    else
        addVal = 0;
    end
    u = u + addVal;
end
v = 0;
for ii=1:length(uniqueVValues)
    n = uniqueVCounts(ii);
    if(k<=n)
        addVal = nchoosek(n,k);
    else
        addVal = 0;
    end
    v = v + addVal;
end

uuCloseToZero = closeToZero(u, len);
vvCloseToZero = closeToZero(v, len);
if( (uuCloseToZero && v>0) || (u>0 && vvCloseToZero) )
    
    % special case of hybrid data
    if(uuCloseToZero)
        continuousRvIndicator = 0;
    else
        continuousRvIndicator = 1;
    end
    numOverlapPtsVec = countOverlaps(U, V, continuousRvIndicator);
    correctionFactor = correctionFactor4(numOverlapPtsVec);
    t = max(u,v)-correctionFactor;
    tau = K/( sqrt(nchoosek(len,2)-t)*sqrt(nchoosek(len,2)-t) );
    
else
    % case of either all continuous or all discrete data
    tau = K/( sqrt(nchoosek(len,2)-u)*sqrt(nchoosek(len,2)-v) );
end

end

function [cf] = correctionFactor4(numOverlapPtsVec)
    meanVal = floor(mean(numOverlapPtsVec));
    if(meanVal<2)
        cf = 0;
    else
        cf = nchoosek(meanVal,2)*length(numOverlapPtsVec);
    end
end

function [out] = closeToZero(in, len)
out = 1;
thresh = 0.02;      % if we are > 2% of length in terms of combinations;
lenFloor = floor(len*thresh);

if(lenFloor>=2)
    cmpVal = nchoosek(lenFloor,2);
else
    cmpVal = 0;
end

if(in>cmpVal)
    out = 0;
end

end

function [numOverlapPtsVec] = countOverlaps(U, V, continuousRvIndicator)
% this function is only called for hybrid data, attempts to correct for
% overestimation of the number of ties in hybrid data

M = length(U);

if(continuousRvIndicator==0)
    % U is the continuous RV
    continuousOutcomes = U;
    discreteOutcomes = V;
    % get the number of unique discrete outcomes
    uniqueDiscreteOutcomes = unique(V);
else
    % V is the continuous RV
    continuousOutcomes = V;
    discreteOutcomes = U;
    % get the number of unique discrete outcomes
    uniqueDiscreteOutcomes = unique(U);
end
uniqueDiscreteOutcomes = sort(uniqueDiscreteOutcomes);  % probably already comes sorted?

% for each unique outcome .. count the overlapping elements.
numOverlapPtsVec = zeros(1,length(uniqueDiscreteOutcomes)-1);
for discreteOutcomesIdx=1:length(uniqueDiscreteOutcomes)-1
    % find the min/max values of the continuous values for this idx and the
    % next
    I = discreteOutcomes==uniqueDiscreteOutcomes(discreteOutcomesIdx);
    relevantContinuousOutcomes_curIdx = continuousOutcomes(I);
    I = discreteOutcomes==uniqueDiscreteOutcomes(discreteOutcomesIdx+1);
    relevantContinuousOutcomes_nextIdx = continuousOutcomes(I);
    
    % compute the number of points which are overlapping
    minCur = min(relevantContinuousOutcomes_curIdx);
    maxCur = max(relevantContinuousOutcomes_curIdx);
    
    numOverlapPoints = length(find(relevantContinuousOutcomes_nextIdx>=minCur & ...
                                   relevantContinuousOutcomes_nextIdx<=maxCur));
                               
    numOverlapPtsVec(discreteOutcomesIdx) = numOverlapPoints/length(relevantContinuousOutcomes_nextIdx)*(M/length(uniqueDiscreteOutcomes));

end

end
