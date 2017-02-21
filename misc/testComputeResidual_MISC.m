clear;
clc;

depType = 'quadratic';

if(strcmpi(depType, 'quadratic'))
    fname = 'quadratic.mat';
elseif(strcmpi(depType, 'circular'))
    fname = 'circular.mat';
elseif(strcmpi(depType, 'sinusoidal'))
    fname = 'sinusoidal.mat';
elseif(strcmpi(depType, 'fourthroot'))
    fname = 'fourthroot.mat';
elseif(strcmpi(depType, 'linear'))
    fname = 'linear.mat';
end

if(ispc)
    % TODO:
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/residualTesting/%s',fname));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/residualTesting/%s',fname));
end

[resid, residAssocIdxs] = ...
    computeResidual_MISC(rectangleCell, ax1pts, ax2pts);