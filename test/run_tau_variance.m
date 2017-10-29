%% Script whch explores how much tau varies as we scan through a pattern
% for different noise levels

clear;
clc;
close all;
dbstop if error;

% find the value of Tau which corresponds to a linear dependency + a
% certain noise level

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec=num_noise_test_min:num_noise_test_max;

MVec = [10 25 50 75 100 250 500 750 1000];
numMCSim = 1000;

resultsVec = zeros(length(MVec),length(noiseVec),numMCSim);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for MIdx=1:length(MVec)
    M = MVec(MIdx);
    for lIdx=1:length(noiseVec)
        l = noiseVec(lIdx);
        dispstat(sprintf('Computing for M=%d noise=%d',M, l));
        tauValTotal = 0;
        parfor mcSimNum=1:numMCSim
            x = rand(M,1);
            y = x + noise*(l/num_noise)*randn(M,1);

            % measure the tau
            tauVal = corr(x,y,'type','kendall');
            resultsVec(MIdx,lIdx,mcSimNum) = tauVal;
        end
    end
end

% save off the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\tau_variance.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/clustering/tau_variance.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/clustering/tau_variance.mat');
end

%% Understand the variance of tau from teh results of the above cell's simulation
clear;
clc;
if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\tau_variance.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/clustering/tau_variance.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/clustering/tau_variance.mat');
end

% MVecToPlot = [10 25 50 75 100 250 500 750 1000];
MVecToPlot = [50 250 500 1000];

f1 = figure(1); 
f2 = figure(2);
f3 = figure(3);
p = numSubplots(length(MVecToPlot));
noiseVecToPlot = noiseVec/10;
labelCell = cell(1,length(noiseVecToPlot));
for ii=1:length(noiseVecToPlot)
    if(ii==1 || ii==11 || ii==21 || ii==31)
        labelCell{ii} = sprintf('%0.01f',noiseVecToPlot(ii));
    else
        labelCell{ii} = '';
    end
end
fontSizeVal = 20;
for ii=1:length(MVecToPlot)
    M = MVecToPlot(ii);
    MIdx = find(M==MVec);

    data = squeeze(resultsVec(MIdx,:,:));
    data = data';  % maek each column the noise
    
%     figure(f1);
%     h = subplot(p(1),p(2),ii);
%     boxplot(data,'Labels',labelCell); 
%     grid on; 
%     title(sprintf('M=%d',M),'FontSize',fontSizeVal);
%     h.FontSize = fontSizeVal;
    
    figure(f2);
    varDataEmpirical = var(data); %compute variance from mean
    meanDataEmpirical = mean(data);
    h = subplot(p(1),p(2),ii);
    varTheoretical1 = ones(1,length(noiseVec))*(2*(2*M+5))./(9*M.*(M-1));
    varTheoretical2 = 4*(M-1.5)/(M*(M-1))*(1-meanDataEmpirical.^2);
    varTheoretical3 = 4*(1-meanDataEmpirical.^2);
    
    stdTheoretical1 = (1-meanDataEmpirical).*sqrt(varTheoretical1);
    stdTheoretical2 = sqrt(varTheoretical2);
    stdTheoretical3 = sqrt(varTheoretical3)/sqrt(M);
    alpha = 0.05;
    
    plot(noiseVecToPlot,sqrt(varDataEmpirical),...
         noiseVecToPlot,stdTheoretical1,...
         noiseVecToPlot,stdTheoretical2,...
         noiseVecToPlot,stdTheoretical3,'o',...
         'LineWidth',2.5);
    grid on; 
    title(sprintf('M=%d',M),'FontSize',fontSizeVal);
    h.FontSize = fontSizeVal;
    
    figure(f3);
    h = subplot(p(1),p(2),ii);
    plot(noiseVecToPlot,stdTheoretical1*4,...
         noiseVecToPlot,stdTheoretical3*norminv(1-alpha),'o',...
         'LineWidth',2.5);
    grid on; 
    title(sprintf('M=%d',M),'FontSize',fontSizeVal);
    h.FontSize = fontSizeVal;
    
end

% figure(f1);
% [~,h] = suplabel('Noise', 'x');
% h.FontSize = fontSizeVal;
% [~,h] = suplabel('\tau','y');
% h.FontSize = fontSizeVal;

figure(f2);
[~,h] = suplabel('Noise', 'x');
h.FontSize = fontSizeVal;
[~,h] = suplabel('\sigma(\tau)','y');
h.FontSize = fontSizeVal;
legend('Empirical','Predicted-1','Predicted-2', 'Predicted-3');

figure(f3);
[~,h] = suplabel('Noise', 'x');
h.FontSize = fontSizeVal;
[~,h] = suplabel('\sigma(\tau)','y');
h.FontSize = fontSizeVal;
legend('Predicted-1','Predicted-3');

%% Check the values of srho, rho, and ktau for linear relationships to see if a correpsondence exists?
clear;
clc;

numMCSim = 500;
M = 500;
lVec = [1,2,3];

resMat = zeros(length(lVec),numMCSim,3);
for ii=1:length(lVec)
    lVal = lVec(ii);
    for mcSimNum=1:numMCSim
        xx = rand(M,1);
        yy = xx + randn(M,1)*lVal;
        rhoVal = corr(xx,yy);
        srhoVal = corr(xx,yy,'type','spearman');
        tauVal = corr(xx,yy,'type','kendall');
        
        resMat(ii,mcSimNum,1) = rhoVal;
        resMat(ii,mcSimNum,2) = srhoVal;
        resMat(ii,mcSimNum,3) = tauVal;
    end
end

% save off the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\linear_relationship_results.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/clustering/linear_relationship_results.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/clustering/linear_relationship_results.mat');
end
