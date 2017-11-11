% script to merge all the power-curves into one
clear;
clc;

load('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
powerTensor_alienware = powerTensor; MVec_alienware = MVec; MIdx_alienware = MIdx;
load('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all_oldwindows.mat');
powerTensor_oldwindows = powerTensor; MVec_oldwindows = MVec; MIdx_oldwindows = MIdx;
load('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all_linux.mat');
powerTensor_linux = powerTensor; MVec_linux = MVec; MIdx_linux = MIdx;

powerTensor_final = powerTensor_alienware;
powerTensor_final(62:62+9-1,:,:,:) = powerTensor_oldwindows(5:13,:,:,:);
powerTensor_final(71:80,:,:,:) = powerTensor_linux(10:-1:1,:,:,:);
powerTensor = powerTensor_final;
MVec = MVec_alienware;
save('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');