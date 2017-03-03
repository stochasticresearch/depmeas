%% Test the Mutual Information based CD metric interface
clear;
clc;

% load the latest version of mutualInfoCD
clear classes;
lnnmod = py.importlib.import_module('lnn');
py.reload(lnnmod);

M = 1000;

x = rand(1,M);
y = rand(1,M);
z = rand(1,M);

lnnCDMetric = lnn_CD(x,y,z);

%% Test Partial Distance Correlation Interface
clear;
clc;

M = 1000;

x = rand(1,M);
y = rand(1,M);
z = rand(1,M);

[pdcor_val] = pdcorr_R(x,y,z)
[pdcor_val, pval] = pdcorr_R(x,y,z)

%% Test Projection Distance Correlation Interface
clear;
clc;

M = 1000;

x = rand(1,M);
y = rand(1,M);
z = rand(1,M);

projcor_val = projcorr_R(x,y,z)