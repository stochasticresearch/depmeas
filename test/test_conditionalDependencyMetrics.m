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

% Independent
x = rand(M,1); y = rand(M,1); z = rand(M,1);
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> Independence\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% Correlated-1
xyz = mvnrnd([0 0 0], [1 0 .3; 0 1 -.1; .3 -.1 1], M);
x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> Correlated-1\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% Correlated-2
xyz = mvnrnd([0 0 0], [1 0 .2; 0 1 -.8; .2 -.8 1], M);
x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> Correlated-2\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% ^structure-1
x = rand(M,1);
gamma = 0.2; eps = randn(M,1);
y = gamma*x + (1-gamma)*eps;
z = gamma*x + (1-gamma)*eps;
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> Structure-1\n');
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
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> Structure-2\n');
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
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> Structure-3\n');
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
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> V-Structure-1\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% vstructure-2
y = rand(M,1); z = rand(M,1);
gamma = 0.2; eps = randn(M,1);
x = gamma*(y+z) + (1-gamma)*eps;
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> V-Structure-2\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);

% vstructure-3
y = rand(M,1); z = rand(M,1);
gamma = 0.2; eps = randn(M,1);
x = gamma*(y+z) + (1-gamma)*eps;
cmVal = cm(x,y,z);
hdVal = hd(x,y,z);
[~, pdcorVal] = pdcov(x,y,z, numR);
hsicVal = hsic(x,y,z,numR);
data = struct(); data.X = x; data.Y = y; data.Z = z; cmiVal = cmi(data);
cmaVal = cassor(data);
rscdmVal = rscdm(x,y,z);
fprintf('>>> V-Structure-3\n');
fprintf('\t cm=%0.03f\n', cmVal);
fprintf('\t hd=%0.03f\n', hdVal);
fprintf('\t pdcor=%0.03f\n', pdcorVal);
fprintf('\t hsic=%0.03f\n', hsicVal);
fprintf('\t cmi=%0.03f\n', cmiVal);
fprintf('\t cma=%0.03f\n', cmaVal);
fprintf('\t rscdm=%0.03f\n', rscdmVal);
