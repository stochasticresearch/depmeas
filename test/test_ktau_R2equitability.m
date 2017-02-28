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

clear;
clc;

rng(123);

% Generate the equitability plot for Kendall's Tau
xi = 0.01:0.01:0.99;

nrow = 99; ncol = 100;

res_f1 = zeros(nrow, ncol);
res_f2 = zeros(nrow, ncol);

R = 0.01;
Rinv = 1/R;
M = 500;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the equitability simulation...\n'),'keepthis','timestamp');

for ii=1:nrow
    dispstat(sprintf('%0.02f %% Complete',ii/nrow*100),'timestamp');
    parfor jj=1:ncol       
        x = rand(M,1)*10+2;
        y = x;
        noise = randn(M,1)*sqrt(var(y)*(Rinv-1));
        y = y+noise;
        res_f1(ii,jj) = rsdm(x,y);
        
        x = rand(M,1)*10+2;
        y = exp(x);
        noise = randn(M,1)*sqrt(var(y)*(Rinv-1));
        y = y+noise;
        res_f2(ii,jj) = rsdm(x,y);
    end
    
    R = R + 0.01;
    Rinv = 1/R;
end

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\equitability_ktau.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/equitability_ktau.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/equitability_ktau');
end

%% load the data and plot
if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\equitability_ktau.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/independence/equitability_ktau.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/independence/equitability_ktau');
end

fontSize = 20;

hold on;
for jj=1:99
    plot(xi, res_f1(:,jj), 'g*');
    plot(xi, res_f2(:,jj), 'b*');
end

title('Equitability Curves for \tau', 'FontSize', fontSize);
grid on;
xlabel('\Phi', 'FontSize', fontSize); ylabel('\tau', 'FontSize', fontSize);
h = legend('Y = X','Y = e^X', 'Location', 'NorthWest');
set(h,'FontSize', fontSize);

%% Plot Figure 5 - which shows the effect of statistical ranks and noise

clear;
clc;

rng(123);

fontSize = 20;

M = 500;
x = rand(M,1)*10+2;
y1 = x; y2 = exp(x);
sigma = 2;
noise = sigma*randn(M,1);

u = pobs(x);
v1 = pobs(y1); v1_n = pobs(y1+noise);
v2 = pobs(y2); v2_n = pobs(y2+noise);

subplot(2,2,1); scatter(x,y1); grid on;
xlabel({'x','(a)'}, 'FontSize', 20); ylabel('y', 'FontSize', 20);
subplot(2,2,2); scatter(u,v1_n); grid on;
xlabel({'u','(b)'}, 'FontSize', 20); ylabel('v', 'FontSize', 20);

subplot(2,2,3); scatter(x,y2); grid on;
xlabel({'x','(c)'}, 'FontSize', 20); ylabel('y', 'FontSize', 20);
subplot(2,2,4); scatter(u,v2_n); grid on;
xlabel({'u','(d)'}, 'FontSize', 20); ylabel('v', 'FontSize', 20);