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
clear;
clc;

M = 500;
numMCSims = 500;
% Test Lower Tail Computation
alphaVec = 0:1:10;
for alpha=alphaVec
    lt_dep_accum = 0;
    for jj=1:numMCSims
        U1 = copularnd('Clayton', alpha, M);
        lt_dep = ltdepest(U1);
        lt_dep_accum = lt_dep_accum + lt_dep;
    end
    lt_dep = lt_dep_accum/numMCSims;
    lt_theo = 2^(-1/alpha);
    fprintf('LT: est=%0.02f actual=%0.02f\n', lt_dep, lt_theo);
%     scatter(U1(:,1),U1(:,2)); grid on;
%     pause;
end

% Test Upper Tail computation
alphaVec = 1:10;
for alpha=alphaVec
    ut_dep_accum = 0;
    for jj=1:numMCSims
        U1 = copularnd('Clayton', alpha, M);
        ut_dep = utdepest(U1);
        ut_dep_accum = ut_dep_accum + ut_dep;
    end
    ut_dep = ut_dep_accum/numMCSims;
    ut_theo = 2-2^(1/alpha);
    fprintf('UT: est=%0.02f actual=%0.02f\n', ut_dep, ut_theo);
%     scatter(U1(:,1),U1(:,2)); grid on;
%     pause;
end
