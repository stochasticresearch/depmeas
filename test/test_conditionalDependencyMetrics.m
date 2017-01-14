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

%% Sanity check of all the measures of conditional association

clear;
clc;

rng(123);
M = 500;
numR = 199;

% ^structure-1
x = rand(M,1);
gamma = 0.2; eps = randn(M,1);
y = gamma*x + (1-gamma)*eps;
z = gamma*x + (1-gamma)*eps;
cmVal = cm(y,z,x);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(y,z,x,numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(y,z,x);
fprintf('>>> ^ Structure-1\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% ^structure-2
x = rand(M,1);
gamma = 0.5; eps = randn(M,1);
y = gamma*x + (1-gamma)*eps;
z = gamma*x + (1-gamma)*eps;
cmVal = cm(y,z,x);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(y,z,x,numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(y,z,x);
fprintf('>>> ^ Structure-2\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% ^structure-3
x = rand(M,1);
gamma = 0.5; eps = randn(M,1);
y = gamma*x + (1-gamma)*eps;
z = gamma*x + (1-gamma)*eps;
cmVal = cm(y,z,x);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(y,z,x,numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(y,z,x);
fprintf('>>> ^ Structure-3\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% vstructure-1
y = rand(M,1); z = rand(M,1);
gamma = 0.2; eps = randn(M,1);
x = gamma*(y+z) + (1-gamma)*eps;
cmVal = cm(y,z,x);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(y,z,x,numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(y,z,x);
fprintf('>>> v Structure-1\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% vstructure-2
y = rand(M,1); z = rand(M,1);
gamma = 0.5; eps = randn(M,1);
x = gamma*(y+z) + (1-gamma)*eps;
cmVal = cm(y,z,x);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(y,z,x,numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(y,z,x);
fprintf('>>> v Structure-2\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% vstructure-3
y = rand(M,1); z = rand(M,1);
gamma = 0.9; eps = randn(M,1);
x = gamma*(y+z) + (1-gamma)*eps;
cmVal = cm(y,z,x);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(y,z,x,numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(y,z,x);
fprintf('>>> v Structure-3\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);
