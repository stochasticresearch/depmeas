function [metric,regionRectangle] = cim(x, y, varargin)
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

if(nargin>=3)
    msi = varargin{1};
else
    msi = 0.015625;
end
[metric,regionRectangle] = cim_cc_mex(x,y,msi);