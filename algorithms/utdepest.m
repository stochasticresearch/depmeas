function [srho_ut] = utdepest(X)
%LTDEPEST - attempts to estimate the upper-tail dependence coefficient in
%terms of Spearman's Rho
% Inputs:
%  X - a matrix of dimensionality [N x D], where N is the # of samples and
%      D is the dimensionality
% Outputs:
%  srho_lt - an estimate of the upper tail dependence coefficient
%
% Ported from ITE Toolbox: https://bitbucket.org/szzoli/ite
%
%**************************************************************************
%* 
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>
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

U = pobs(X);
d = size(U,2);
p = 0.5;

c = mean(prod(1 - max(U',1-p)));
c1 = (p*(2-p)/2)^d;
c2 = p^d * (d+1-p*d) / (d+1);

srho_ut = (c-c1) / (c2-c1);

end