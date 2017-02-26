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

%% Show w/ a simple plot that Kendall's Tau is not equitable

noisePowerVec = 0:.1:3;
M = 500;
numDepTypes = 2;    % 1 = linear, 2 = exponential
tauResMatrix = zeros(numDepTypes, length(noisePowerVec));

for noisePowerIdx=1:length(noisePowerVec)
    noisePower = noisePowerVec(noisePowerIdx);
    x = rand(M,1)*10+2;
    noise = randn(M,1)*noisePower;
    
    y1 = x + noise;
    y2 = exp(x) + noise;
    
    tau1 = corr(x,y1,'type','kendall');
    tau2 = corr(x,y2,'type','kendall');
    
    tauResMatrix(1,noisePowerIdx) = tau1;
    tauResMatrix(2,noisePowerIdx) = tau2;
end

plot(noisePowerVec, tauResMatrix);

%% Generate Figure 5
clear;
clc;

M = 1000;
x = rand(M,1)*8+2;
y1 = x;
y2 = exp(x);
noise = randn(M,1)*2;

u = pobs(x);
v1 = pobs(y1+noise);
v2 = pobs(y2+noise);

subplot(2,2,1);
scatter(x,y1); grid on; xlabel({'x', '(a)'}, 'FontSize', 20); ylabel('y', 'FontSize', 20);
subplot(2,2,2);
scatter(u,v1); grid on; xlabel({'u', '(b)'}, 'FontSize', 20); ylabel('v', 'FontSize', 20);
subplot(2,2,3);
scatter(x,y2); grid on; xlabel({'x', '(c)'}, 'FontSize', 20); ylabel('y', 'FontSize', 20);
subplot(2,2,4);
scatter(u,v2); grid on; xlabel({'u', '(d)'}, 'FontSize', 20); ylabel('v', 'FontSize', 20);
