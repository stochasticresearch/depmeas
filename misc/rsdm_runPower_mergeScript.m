%% Merge the CoS/cCorr/TICe data
clear;
clc;
dbstop if error;

loRangeFilename = 'rsdmPower_CoS_cCorr_ticE_M_25_825.mat';
hiRangeFilename = 'rsdmPower_CoS_cCorr_ticE_M_850_1500.mat';
outputFilename  = 'rsdmPower_CoS_cCorr_ticE_M_25_1500.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, loRangeFilename);
f2 = fullfile(folder, hiRangeFilename);

load(f2);
ccorrPower_toMerge = ccorrPower;
cosPower_toMerge = cosPower;
ticePower_toMerge = ticePower;
load(f1);
mergeStartIdx = length(25:25:825);
numToMerge = length(850:25:1500);

% merge the powers
for typ=1:numDepTests
    for l=num_noise_test_min:num_noise_test_max
        for mm=1:numToMerge
            % for the submerge
            ccorrPower(typ,l,mm+mergeStartIdx) = ccorrPower_toMerge(typ,l,mm+mergeStartIdx);
            cosPower(typ,l,mm+mergeStartIdx) = cosPower_toMerge(typ,l,mm+mergeStartIdx);
            ticePower(typ,l,mm+mergeStartIdx) = ticePower_toMerge(typ,l,mm+mergeStartIdx);
        end
    end
end

% save as rsdmPower.mat
finalOutputFile = fullfile(folder, outputFilename);
clearvars loRangeFilename hiRangeFilename outputFilename
clearvars folder loFilename hiFilename
clearvars ccorrPower_toMerge cosPower_toMerge ticePower_toMerge
clearvars mergeStartIdx numToMerge 
save(finalOutputFile);
%% RSDM Power Merge Script
clear;
clc;
dbstop if error;

% loRangeFilename = 'rsdmPower_M_25_400.mat';
% hiRangeFilename = 'rsdmPower_M_425_750.mat';
% outputFilename  = 'rsdmPower_M_25_750.mat';
loRangeFilename = 'rsdmPower_M_25_750.mat';
hiRangeFilename = 'rsdmPower_M_775_1500.mat';
outputFilename  = 'rsdmPower.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, loRangeFilename);
f2 = fullfile(folder, hiRangeFilename);

load(f2);
corrPower_toMerge = corrPower;
dcorrPower_toMerge = dcorrPower;
micePower_toMerge = micePower;
rdcPower_toMerge = rdcPower;
rsdmPower_toMerge = rsdmPower;
load(f1);
% mergeStartIdx = length(25:25:400);
% numToMerge = length(425:25:750);
% M_vec = 25:25:750;
mergeStartIdx = length(25:25:750);
numToMerge = length(775:25:1500);
M_vec = 25:25:1500;


% Merge the M_vec
% no need to merge this, b/c we controlled the runs by teh idxs directly

% merge the powers
for typ=1:numDepTests
    for l=num_noise_test_min:num_noise_test_max
        for mm=1:numToMerge
            % for the submerge
            corrPower(typ,l,mm+mergeStartIdx) = corrPower_toMerge(typ,l,mm+mergeStartIdx);
            dcorrPower(typ,l,mm+mergeStartIdx) = dcorrPower_toMerge(typ,l,mm+mergeStartIdx);
            micePower(typ,l,mm+mergeStartIdx) = micePower_toMerge(typ,l,mm+mergeStartIdx);
            rdcPower(typ,l,mm+mergeStartIdx) = rdcPower_toMerge(typ,l,mm+mergeStartIdx);
            rsdmPower(typ,l,mm+mergeStartIdx) = rsdmPower_toMerge(typ,l,mm+mergeStartIdx);
        end
    end
end

% save as rsdmPower.mat
finalOutputFile = fullfile(folder, outputFilename);
clearvars loRangeFilename hiRangeFilename outputFilename
clearvars folder loFilename hiFilename
clearvars corrPower_toMerge dcorrPower_toMerge micePower_toMerge rdcPower_toMerge rsdmPower_toMerge
clearvars mergeStartIdx numToMerge 
save(finalOutputFile);

%% Merge the powers of rsdm/dCorr/MICe/rdc/corr with CoS/cCorr/TICe

clear;
clc;
dbstop if error;

f1name = 'power_rsdm_dCorr_rdc_corr_mice.mat';
f2name = 'power_CoS_cCorr_ticE.mat';
outputFilename = 'power_all.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f2);
load(f1);

% save as rsdmPower.mat
finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);

%% Merge cosdv power into power_all

clear;
clc;
dbstop if error;

f1name = 'CoSPower_M_25_1500.mat';
f2name = 'power_all.mat';
outputFilename = 'power_all.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f1);
cosdvPower = cosPower;
load(f2);
cosfpower = cosPower;       % save off old simulations for cosf
cosPower = cosdvPower;      % overwrite old cosf w/ cosdv simulations

finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);

%% Merge the cosdv power into power_M_500 (i.e. extract the relevant vector)

clear;
clc;
dbstop if error;

f1name = 'power_all.mat';
f2name = 'power_M_500.mat';
outputFilename = 'power_M_500.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f1);
clearvars -except cosPower f1name f2name outputFilename folder f1 f2
cosdvPower = cosPower;
load(f2);
cosfpower = cosPower;   % save off old simulations for cosf
cosPower = squeeze(cosdvPower(:,:,20));  % overwrite old cosf w/ cosdv simulations, grab M=500 data
clearvars cosdvPower    % no need to store redundant information, makes for confusion later :D

finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);

%% Merge knn power into power_M_500

clear;
clc;
dbstop if error;

f1name = 'rsdmPower_knn_M_500.mat';
f2name = 'power_M_500.mat';
outputFilename = 'power_M_500.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f1);
clearvars -except knn1Power knn6Power knn20Power f1name f2name outputFilename folder f1 f2
load(f2);

finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);

%% Merge the ite power into power_M_500

clear;
clc;
dbstop if error;

f1name = 'rsdmPower_ite_M_500.mat';
f2name = 'power_M_500.mat';
outputFilename = 'power_M_500.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f1);
clearvars -except shapPower shvmePower f1name f2name outputFilename folder f1 f2
load(f2);

finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);

%% Merge the RSDM only runs

clear;
clc;
dbstop if error;

f1name = 'rsdmPower_rsdm_M_25_650_overfitstrategy1.mat';
f2name = 'rsdmPower_rsdm_M_675_1050_overfitstrategy1.mat';
outputFilename = 'rsdmPower_rsdm_25_1050_overfitstrategy1.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f2);
rsdmPower_toMerge = rsdmPower;
load(f1);

mergeStartIdx = length(25:25:650);
numToMerge = length(675:25:1050);
M_vec = 25:25:1500;

% merge the powers
for typ=1:numDepTests
    for l=num_noise_test_min:num_noise_test_max
        for mm=1:numToMerge
            % for the submerge
            rsdmPower(typ,l,mm+mergeStartIdx) = rsdmPower_toMerge(typ,l,mm+mergeStartIdx);
        end
    end
end

% save the output
finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
clearvars rsdmPower_toMerge
clearvars mergeStartIdx numToMerge 
save(finalOutputFile);

%% RSDM -- merge the backwards runs into the forward runs now

clear;
clc;
dbstop if error;

f1name = 'rsdmPower_rsdm_25_1050_overfitstrategy1.mat';
f2name = 'rsdmPower_rsdm_M_1500_1075_overfitstrategy1.mat';
outputFilename = 'rsdmPower_rsdm_M_25_1500_overfitstrategy1.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f2);
rsdmPower_toMerge = rsdmPower;
load(f1);

mergeStartIdx = length(25:25:1050);
numToMerge = length(1075:25:1500);
M_vec = 25:25:1500;

% merge the powers
for typ=1:numDepTests
    for l=num_noise_test_min:num_noise_test_max
        rsdmPowerIdx = length(M_vec);
        for mm=1:numToMerge
            % for the submerge
            fprintf('rsdmPowerIdx=%d mm=%d\n', rsdmPowerIdx, mm);
            rsdmPower(typ,l,rsdmPowerIdx) = rsdmPower_toMerge(typ,l,mm);
            rsdmPowerIdx = rsdmPowerIdx - 1;    % b/c we are merging backwards runs
        end
    end
end

% save the output
finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
clearvars rsdmPower_toMerge
clearvars mergeStartIdx numToMerge 
save(finalOutputFile);

%% Sanity Check of whether the overfit strategys are helping or not ...

clear;
clc;

f1name = 'rsdmPower_rsdm_M_25_1500_overfitstrategy1.mat';
f2name = 'power_all.mat';

if(ispc)
    folder1 = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence\\backup';
    folder2 = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder1 = '/Users/Kiran/ownCloud/PhD/sim_results/independence/backup';
    folder2 = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder1 = '/home/kiran/ownCloud/PhD/sim_results/independence/backup';
    folder2 = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder1, f1name);
f2 = fullfile(folder2, f2name);

load(f1);
rsdmPower_overfitStrategy1 = rsdmPower;
load(f2);

% first do a sanity check to see if the powers are significantly different,
% if not, we wont bother with the merge ...
% process these results to search for the minimum M that gives us the
% threshold power
rsdmPower = squeeze(rsdmPower(:,:,20));
rsdmOverfitStrategy1Power = squeeze(rsdmPower_overfitStrategy1(:,:,20));

% inlet plot configuration
M_inlet = 200; M = 500;
if(M==500)
    inset_bufX = 0.0005; inset_bufY = 0.002;
else
    inset_bufX = 0.15; inset_bufY = 0.26;
end

inset_width = 0.1; inset_height = 0.08;

noiseVec = (num_noise_test_min:num_noise_test_max)/10;
figure;
h1 = subplot(2,2,1);
hh1 = plot(noiseVec, rsdmPower(1,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(1,num_noise_test_min:num_noise_test_max), '+-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h1.FontSize = 20; 
loc_inset = [h1.Position(1)+inset_bufX h1.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 

h2 = subplot(2,2,2);
hh2 = plot(noiseVec, rsdmPower(2,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(2,num_noise_test_min:num_noise_test_max), '+-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h2.FontSize = 20; 
loc_inset = [h2.Position(1)+inset_bufX h2.Position(2)+inset_bufY inset_width inset_height];
ax2 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 4*(tmp1-.5).^2;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax2.Box = 'on'; ax2.XTick = []; ax2.YTick = [];
ax2.XLim = [min(tmp1) max(tmp1)];
ax2.YLim = [min(tmp2) max(tmp2)];
hh2(1).LineWidth = 1.5; 
hh2(2).LineWidth = 1.5; 

h3 = subplot(2,2,3); 
hh3 = plot(noiseVec, rsdmPower(3,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(3,num_noise_test_min:num_noise_test_max), '+-.');   
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h3.FontSize = 20; 
loc_inset = [h3.Position(1)+inset_bufX h3.Position(2)+inset_bufY inset_width inset_height];
ax3 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 128*(tmp1-1/3).^3-48*(tmp1-1/3).^3-12*(tmp1-1/3);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax3.Box = 'on'; ax3.XTick = []; ax3.YTick = [];
ax3.XLim = [min(tmp1) max(tmp1)];
ax3.YLim = [min(tmp2) max(tmp2)];
hh3(1).LineWidth = 1.5; 
hh3(2).LineWidth = 1.5; 

h4 = subplot(2,2,4); 
hh4 = plot(noiseVec, rsdmPower(4,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(4,num_noise_test_min:num_noise_test_max), '+-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h4.FontSize = 20; 
loc_inset = [h4.Position(1)+inset_bufX h4.Position(2)+inset_bufY inset_width inset_height];
ax4 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(4*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax4.Box = 'on'; ax4.XTick = []; ax4.YTick = [];
ax4.XLim = [min(tmp1) max(tmp1)];
ax4.YLim = [min(tmp2) max(tmp2)];
hh4(1).LineWidth = 1.5; 
hh4(2).LineWidth = 1.5; 

figure;
h5 = subplot(2,2,1); 
hh5 = plot(noiseVec, rsdmPower(5,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(5,num_noise_test_min:num_noise_test_max), '+-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Power'); grid on;
h5.FontSize = 20; 
loc_inset = [h5.Position(1)+inset_bufX h5.Position(2)+inset_bufY inset_width inset_height];
ax5 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(16*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax5.Box = 'on'; ax5.XTick = []; ax5.YTick = [];
ax5.XLim = [min(tmp1) max(tmp1)];
ax5.YLim = [min(tmp2) max(tmp2)];
hh5(1).LineWidth = 1.5; 
hh5(2).LineWidth = 1.5; 

h6 = subplot(2,2,2); 
hh6 = plot(noiseVec, rsdmPower(6,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(6,num_noise_test_min:num_noise_test_max), '+-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Power'); grid on;
h6.FontSize = 20; 
loc_inset = [h6.Position(1)+inset_bufX h6.Position(2)+inset_bufY inset_width inset_height];
ax6 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1.^(1/4);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax6.Box = 'on'; ax6.XTick = []; ax6.YTick = [];
ax6.XLim = [min(tmp1) max(tmp1)];
ax6.YLim = [min(tmp2) max(tmp2)];
hh6(1).LineWidth = 1.5; 
hh6(2).LineWidth = 1.5; 

h7 = subplot(2,2,3); 
hh7 = plot(noiseVec, rsdmPower(7,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(7,num_noise_test_min:num_noise_test_max), '+-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Power'); grid on;
h7.FontSize = 20; 
loc_inset = [h7.Position(1)+inset_bufX h7.Position(2)+inset_bufY inset_width inset_height];
ax7 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet/2);
tmp2 = (sqrt(1 - (2*tmp1 - 1).^2));
tmp3 = -(sqrt(1 - (2*tmp1 - 1).^2));
plot(tmp1,tmp2, 'k', 'LineWidth', 2); hold on;
plot(tmp1,tmp3, 'k', 'LineWidth', 2); 
ax7.Box = 'on'; ax7.XTick = []; ax7.YTick = [];
ax7.XLim = [min(tmp1) max(tmp1)];
ax7.YLim = [min(tmp3) max(tmp2)];
hh7(1).LineWidth = 1.5; 
hh7(2).LineWidth = 1.5; 

h8 = subplot(2,2,4); 
hh8 = plot(noiseVec, rsdmPower(8,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdmOverfitStrategy1Power(8,num_noise_test_min:num_noise_test_max), '+-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
h8.FontSize = 20; 
legend('RSDM', 'RSDM (OS 1)');  % manually move this using the mouse to a
                                                        % good location
xlabel('Noise Level'); ylabel('Power'); grid on;
loc_inset = [h8.Position(1)+inset_bufX h8.Position(2)+inset_bufY inset_width inset_height];
ax8 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = (tmp1 > 0.5);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax8.Box = 'on'; ax8.XTick = []; ax8.YTick = [];
ax8.XLim = [min(tmp1) max(tmp1)];
hh8(1).LineWidth = 1.5; 
hh8(2).LineWidth = 1.5; 

%% Merge the ite power's @ M=100 & RSDM @ M=100

clear;
clc;
dbstop if error;

f1name = 'rsdmPower_MI_M_100.mat';
f2name = 'power_all.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f2);
rsdmPower = squeeze(rsdmPower(:,:,4));  %  grab RSDM M=100 data
clearvars -except rsdmPower f1

load(f1);
clearvars -except rsdmPower knn1Power knn6Power knn20Power shapPower shvmePower

M = 100;
num_noise_test_min = 1;
num_noise_test_max = 30;

outputFilename = 'power_M_100.mat';
if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

finalOutputFile = fullfile(folder, outputFilename);
save(finalOutputFile);