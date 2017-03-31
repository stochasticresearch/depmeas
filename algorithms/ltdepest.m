function [srho_lt] = ltdepest(X)
%LTDEPEST - attempts to estimate the lower-tail dependence coefficient in
%terms of Spearman's Rho
% Inputs:
%  X - a matrix of dimensionality [N x D], where N is the # of samples and
%      D is the dimensionality
% Outputs:
%  srho_lt - an estimate of the lower tail dependence coefficient
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
U = U';

p = 0.5;

c1 = (p^2/2)^d;
c2 = p^(d+1)/(d+1);

srho_lt = (mean(prod(max(p-U,0))) - c1) / (c2 - c1);

end