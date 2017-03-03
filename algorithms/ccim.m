function [ccimVal, RxAligned, RyAligned] = ccim(x, y, z, varargin)
%CCIM - a copula based test of conditional dependence
% Inputs:
%  x - the first random variable
%  y - the second random variable
%  z - the third random variable, which x and y will be conditioned upon
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
%
% Outputs
%  metric - a real-valued number between 0 and 1 which represents the
%           statistical dependence between x|z and y|z
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

% overwrite defaults with user-inputted values
nVarargin = length(varargin);
switch nVarargin
    case 0
        [~, ~, Rx, RxIdxs] = cim(z, x);
        [~, ~, Ry, RyIdxs] = cim(z, y);
    case 1
        minscanincr = varargin{1};
        
        [~, ~, Rx, RxIdxs] = cim(z, x, minscanincr);
        [~, ~, Ry, RyIdxs] = cim(z, y, minscanincr);
    case 2
        minscanincr = varargin{1};
        diffthresh = varargin{2};
        
        [~, ~, Rx, RxIdxs] = cim(z, x, minscanincr, diffthresh);
        [~, ~, Ry, RyIdxs] = cim(z, y, minscanincr, diffthresh);
    otherwise
        % the >=3 case
        minscanincr = varargin{1};
        diffthresh = varargin{2};
        alpha = varargin{3};
        
        [~, ~, Rx, RxIdxs] = cim(z, x, minscanincr, diffthresh, alpha);
        [~, ~, Ry, RyIdxs] = cim(z, y, minscanincr, diffthresh, alpha);
end

% align Rx and Ry
RxStacked = []; RyStacked = [];
RxIdxsStacked = []; RyIdxsStacked = [];
for ii=1:length(RxIdxs)
    tmpIdxs = RxIdxs{ii};
    tmpVals = Rx{ii};
    for jj=1:length(tmpIdxs)
        RxIdxsStacked = [RxIdxsStacked; tmpIdxs{jj}];
        RxStacked = [RxStacked; tmpVals{jj}];
    end
end
for ii=1:length(RyIdxs)
    tmpIdxs = RyIdxs{ii};
    tmpVals = Ry{ii};
    for jj=1:length(tmpIdxs)
        RyIdxsStacked = [RyIdxsStacked; tmpIdxs{jj}];
        RyStacked = [RyStacked; tmpVals{jj}];
    end
end
[~,Ix] = sort(RxIdxsStacked);
[~,Iy] = sort(RyIdxsStacked);

RxAligned = RxStacked(Ix);
RyAligned = RyStacked(Iy);


% TODO: some more processing?
%   1.) ignore box edge points?
%   2.) if magnitude of all points is really small, then we should consider
%       them independent ... is that statistically OK?

switch nVarargin
    case 0
        ccimVal = cim(RxAligned, RyAligned);
    case 1
        minscanincr = varargin{1};
        ccimVal = cim(RxAligned, RyAligned, minscanincr);
    case 2
        minscanincr = varargin{1};
        diffthresh = varargin{2};
        ccimVal = cim(RxAligned, RyAligned, minscanincr, diffthresh);
    otherwise
        % the >=3 case
        minscanincr = varargin{1};
        diffthresh = varargin{2};
        alpha = varargin{3};
        ccimVal = cim(RxAligned, RyAligned, minscanincr, diffthresh, alpha);
end

end