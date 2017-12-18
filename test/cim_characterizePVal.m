%% Characterize null distribution (X indep Y) experimentally for \hat{CIM}
clear;
clc;

% setup and start the parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 5;
saveProfile(myCluster);
p = gcp

rng(1234);

nsim = 1000;
M_vec = 100:100:2500;
initLen = length(M_vec);
M_vec = [M_vec 5000 10000];

xMin = 0; xMax = 1;
yMin = 0; yMax = 1;

minScanIncr = 0.015625;

FIT_PLOTS = 0;

rsdmNullDistributionResultsContinuous = zeros(nsim, length(M_vec));
rsdmNullDistributionResultsDiscrete = zeros(nsim, length(M_vec));
rsdmNullDistributionResultsHybrid1 = zeros(nsim, length(M_vec));
rsdmNullDistributionResultsHybrid2 = zeros(nsim, length(M_vec));

numDiscreteIntervals = 4;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for ii=1:nsim
    dispstat(sprintf('Computing for sim=%d/%d',ii,nsim),'keepthis', 'timestamp');
    parfor jj=1:length(M_vec)
        M = M_vec(jj);
        
        % continuous create independent x & y
        x = rand(M,1)*(xMax-xMin)+xMin;
        y = rand(M,1)*(yMax-yMin)+yMin;
        
        % discrete independent x & y
        x_discrete = discretizeRv(x,numDiscreteIntervals)';
        y_discrete = discretizeRv(y,numDiscreteIntervals)';
    
        % compute CIM
        rsdmNullDistributionResultsContinuous(ii,jj) = cim(x,y,minScanIncr);
        rsdmNullDistributionResultsHybrid1(ii,jj)    = cim(x_discrete,y,minScanIncr);
        rsdmNullDistributionResultsHybrid2(ii,jj)    = cim(x,y_discrete,minScanIncr);
        rsdmNullDistributionResultsDiscrete(ii,jj)   = cim(x_discrete,y_discrete,minScanIncr);
    end
end

% plot distribution of CIM under the null distribution 
M_vec_toPlot = M_vec(1:initLen);
legendCell = cell(1,length(M_vec_toPlot));
subplot(2,2,1);
for ii=1:length(M_vec_toPlot)
    [f,xi] = ksdensity(rsdmNullDistributionResultsContinuous(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec_toPlot(ii));
end
grid on;
legend(legendCell);
title('Distribution of CIM_{approx}, X-C,Y-C');

subplot(2,2,2);
for ii=1:length(M_vec_toPlot)
    [f,xi] = ksdensity(rsdmNullDistributionResultsHybrid1(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec_toPlot(ii));
end
grid on;
legend(legendCell);
title('Distribution of CIM_{approx}, X-D,Y-C');

subplot(2,2,3);
for ii=1:length(M_vec_toPlot)
    [f,xi] = ksdensity(rsdmNullDistributionResultsHybrid2(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec_toPlot(ii));
end
grid on;
legend(legendCell);
title('Distribution of CIM_{approx}, X-C,Y-D');

subplot(2,2,4);
for ii=1:length(M_vec_toPlot)
    [f,xi] = ksdensity(rsdmNullDistributionResultsDiscrete(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec_toPlot(ii));
end
grid on;
legend(legendCell);
title('Distribution of CIM_{approx}, X-D,Y-D');

D_continuous_cell = cell(1,length(M_vec));  PD_continuous_cell = cell(1,length(M_vec));
D_hybrid1_cell = cell(1,length(M_vec));  PD_hybrid1_cell = cell(1,length(M_vec));
D_hybrid2_cell = cell(1,length(M_vec));  PD_hybrid2_cell = cell(1,length(M_vec));
D_discrete_cell = cell(1,length(M_vec));  PD_discrete_cell = cell(1,length(M_vec));
idx = 1;
for ii=1:length(M_vec)
    [D, PD_continuous] = allfitdist(rsdmNullDistributionResultsContinuous(:,ii), 'PDF');
    D_continuous_cell{idx} = D;  PD_continuous_cell{idx} = PD_continuous;
    
    [D, PD_continuous] = allfitdist(rsdmNullDistributionResultsHybrid1(:,ii), 'PDF');
    D_hybrid1_cell{idx} = D;  PD_hybrid1_cell{idx} = PD_continuous;
    
    [D, PD_continuous] = allfitdist(rsdmNullDistributionResultsHybrid2(:,ii), 'PDF');
    D_hybrid2_cell{idx} = D;  PD_hybrid2_cell{idx} = PD_continuous;
    
    [D, PD_continuous] = allfitdist(rsdmNullDistributionResultsDiscrete(:,ii), 'PDF');
    D_discrete_cell{idx} = D;  PD_discrete_cell{idx} = PD_continuous;
    
    idx = idx + 1;
end
if(~FIT_PLOTS)
    close all;      % close the generated plots
end

% for each PD type, compute the total BIC score for all sample sizes, and
% choose the best one in that fashion
distributions = {'Beta', 'Birnbaum-Saunders', 'Exponential', ...
                 'Extreme value', 'Gamma', 'Generalized extreme value', ...
                 'Generalized Pareto', 'Inverse Gaussian', 'Logistic', ...
                 'Log-logistic', 'Lognormal', 'Nakagami', 'Normal', ...
                 'Rayleigh', 'Rician', 't location-scale', 'Weibull'};

distScoresContinuous = zeros(4,length(distributions));
distScoresHybrid1 = zeros(4,length(distributions));
distScoresHybrid2 = zeros(4,length(distributions));
distScoresDiscrete = zeros(4,length(distributions));
for ii=1:length(distributions)
    dist = distributions{ii};
    % find this distribution in the fit and store the BIC, AIC, AICc scores
    % for all M
    NLogL_continuous = 0; BIC_continuous = 0; AIC_continuous = 0; AICc_continuous = 0;
    NLogL_hybrid1 = 0; BIC_hybrid1 = 0; AIC_hybrid1 = 0; AICc_hybrid1 = 0;
    NLogL_hybrid2 = 0; BIC_hybrid2 = 0; AIC_hybrid2 = 0; AICc_hybrid2 = 0;
    NLogL_discrete = 0; BIC_discrete = 0; AIC_discrete = 0; AICc_discrete = 0;
    
    for jj=1:length(M_vec)
        D = D_continuous_cell{jj};
        PD_continuous = PD_continuous_cell{jj};
        % find the distribution
        distFound = 0;
        for kk=1:length(PD_continuous)
            if(strcmpi(PD_continuous{kk}.DistributionName, dist))
                distFound = 1;
                break;
            end
        end
        if(distFound)
            if(isreal(D(kk).NLogL))
                NLogL_continuous = NLogL_continuous + D(kk).NLogL;
            else
                NLogL_continuous = 0;
            end
            if(isreal(D(kk).BIC))
                BIC_continuous = BIC_continuous + D(kk).BIC;
            else
                BIC_continuous = 0;
            end
            if(isreal(D(kk).AIC))
                AIC_continuous = AIC_continuous + D(kk).AIC;
            else
                AIC_continuous = 0;
            end
            if(isreal(D(kk).AICc))
                AICc_continuous = AICc_continuous + D(kk).AICc;
            else
                AICc_continuous = 0;
            end
        else
            NLogL_continuous = 0;
            BIC_continuous = 0;
            AIC_continuous = 0;
            AICc_continuous = 0;
        end
        
        D = D_hybrid1_cell{jj};
        PD_continuous = PD_hybrid1_cell{jj};
        % find the distribution
        distFound = 0;
        for kk=1:length(PD_continuous)
            if(strcmpi(PD_continuous{kk}.DistributionName, dist))
                distFound = 1;
                break;
            end
        end
        if(distFound)
            if(isreal(D(kk).NLogL))
                NLogL_hybrid1 = NLogL_hybrid1 + D(kk).NLogL;
            else
                NLogL_hybrid1 = 0;
            end
            if(isreal(D(kk).BIC))
                BIC_hybrid1 = BIC_hybrid1 + D(kk).BIC;
            else
                BIC_hybrid1 = 0;
            end
            if(isreal(D(kk).AIC))
                AIC_hybrid1 = AIC_hybrid1 + D(kk).AIC;
            else
                AIC_hybrid1 = 0;
            end
            if(isreal(D(kk).AICc))
                AICc_hybrid1 = AICc_hybrid1 + D(kk).AICc;
            else
                AICc_hybrid1 = 0;
            end
        else
            NLogL_hybrid1 = 0;
            BIC_hybrid1 = 0;
            AIC_hybrid1 = 0;
            AICc_hybrid1 = 0;
        end
        
        D = D_hybrid2_cell{jj};
        PD_continuous = PD_hybrid2_cell{jj};
        % find the distribution
        distFound = 0;
        for kk=1:length(PD_continuous)
            if(strcmpi(PD_continuous{kk}.DistributionName, dist))
                distFound = 1;
                break;
            end
        end
        if(distFound)
            if(isreal(D(kk).NLogL))
                NLogL_hybrid2 = NLogL_hybrid2 + D(kk).NLogL;
            else
                NLogL_hybrid2 = 0;
            end
            if(isreal(D(kk).BIC))
                BIC_hybrid2 = BIC_hybrid2 + D(kk).BIC;
            else
                BIC_hybrid2 = 0;
            end
            if(isreal(D(kk).AIC))
                AIC_hybrid2 = AIC_hybrid2 + D(kk).AIC;
            else
                AIC_hybrid2 = 0;
            end
            if(isreal(D(kk).AICc))
                AICc_hybrid2 = AICc_hybrid2 + D(kk).AICc;
            else
                AICc_hybrid2 = 0;
            end
        else
            NLogL_hybrid2 = 0;
            BIC_hybrid2 = 0;
            AIC_hybrid2 = 0;
            AICc_hybrid2 = 0;
        end
        
        D = D_discrete_cell{jj};
        PD_continuous = PD_discrete_cell{jj};
        % find the distribution
        distFound = 0;
        for kk=1:length(PD_continuous)
            if(strcmpi(PD_continuous{kk}.DistributionName, dist))
                distFound = 1;
                break;
            end
        end
        if(distFound)
            if(isreal(D(kk).NLogL))
                NLogL_discrete = NLogL_discrete + D(kk).NLogL;
            else
                NLogL_discrete = 0;
            end
            if(isreal(D(kk).BIC))
                BIC_discrete = BIC_discrete + D(kk).BIC;
            else
                BIC_discrete = 0;
            end
            if(isreal(D(kk).AIC))
                AIC_discrete = AIC_discrete + D(kk).AIC;
            else
                AIC_discrete = 0;
            end
            if(isreal(D(kk).AICc))
                AICc_discrete = AICc_discrete + D(kk).AICc;
            else
                AICc_discrete = 0;
            end
        else
            NLogL_discrete = 0;
            BIC_discrete = 0;
            AIC_discrete = 0;
            AICc_discrete = 0;
        end
    end
    
    distScoresContinuous(1,ii) = NLogL_continuous;
    distScoresContinuous(2,ii) = BIC_continuous;
    distScoresContinuous(3,ii) = AIC_continuous;
    distScoresContinuous(4,ii) = AICc_continuous;
    
    distScoresHybrid1(1,ii) = NLogL_hybrid1;
    distScoresHybrid1(2,ii) = BIC_hybrid1;
    distScoresHybrid1(3,ii) = AIC_hybrid1;
    distScoresHybrid1(4,ii) = AICc_hybrid1;
    
    distScoresHybrid2(1,ii) = NLogL_hybrid2;
    distScoresHybrid2(2,ii) = BIC_hybrid2;
    distScoresHybrid2(3,ii) = AIC_hybrid2;
    distScoresHybrid2(4,ii) = AICc_hybrid2;
    
    distScoresDiscrete(1,ii) = NLogL_discrete;
    distScoresDiscrete(2,ii) = BIC_discrete;
    distScoresDiscrete(3,ii) = AIC_discrete;
    distScoresDiscrete(4,ii) = AICc_discrete;
end

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_nullDistribution.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_nullDistribution.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/cim_nullDistribution.mat');
end

fprintf('*************** X & Y CONTINUOUS ****************\n');
% Sort by NLogL
[~,I] = sort(distScoresContinuous(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresContinuous(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresContinuous(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresContinuous(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\n');

fprintf('*************** X & Y HYBRID 1 ****************\n');
% Sort by NLogL
[~,I] = sort(distScoresHybrid1(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresHybrid1(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresHybrid1(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresHybrid1(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\n');

fprintf('*************** X & Y HYBRID 2 ****************\n');
% Sort by NLogL
[~,I] = sort(distScoresHybrid2(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresHybrid2(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresHybrid2(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresHybrid2(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\n');

fprintf('*************** X & Y DISCRETE ****************\n');
% Sort by NLogL
[~,I] = sort(distScoresDiscrete(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresDiscrete(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresDiscrete(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresDiscrete(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\n');

%%
clear;
clc;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_nullDistribution.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_nullDistribution.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/independence/cim_nullDistribution.mat');
end


% From the above analysis, the Beta distribution seems 
% to fit best ... Q-Q Plots

% QQ Plot w/ best fit for M=100 and M=1000
pdObjsContinuous = cell(1,length(M_vec));
alphaVecContinuous = zeros(1,length(M_vec));
alphaVecHybrid1 = zeros(1,length(M_vec));
alphaVecHybrid2 = zeros(1,length(M_vec));
alphaVecDiscrete = zeros(1,length(M_vec));
betaVecContinuous = zeros(1,length(M_vec));
betaVecHybrid1 = zeros(1,length(M_vec));
betaVecHybrid2 = zeros(1,length(M_vec));
betaVecDiscrete = zeros(1,length(M_vec));
for ii=1:length(M_vec)
    M = M_vec(ii);
    % look for the Inverse Gaussian Distribution in the correct cell array
    PD_continuous = PD_continuous_cell(ii); PD_continuous = PD_continuous{1};
    PD_hybrid1 = PD_hybrid1_cell(ii); PD_hybrid1 = PD_hybrid1{1};
    PD_hybrid2 = PD_hybrid2_cell(ii); PD_hybrid2 = PD_hybrid2{1};
    PD_discrete = PD_discrete_cell(ii); PD_discrete = PD_discrete{1};
    for jj=1:length(PD_continuous)
        if(strcmpi('Beta', PD_continuous{jj}.DistributionName))
            pdContinuous = PD_continuous{jj};
            break;
        end
    end
    for jj=1:length(PD_hybrid1)
        if(strcmpi('Beta', PD_hybrid1{jj}.DistributionName))
            pdHybrid1 = PD_hybrid1{jj};
            break;
        end
    end
    for jj=1:length(PD_hybrid2)
        if(strcmpi('Beta', PD_hybrid2{jj}.DistributionName))
            pdHybrid2 = PD_hybrid2{jj};
            break;
        end
    end
    for jj=1:length(PD_discrete)
        if(strcmpi('Beta', PD_discrete{jj}.DistributionName))
            pdDiscrete = PD_discrete{jj};
            break;
        end
    end
    
    pdObjsContinuous{ii} = pdContinuous;
    alphaVecContinuous(ii) = pdContinuous.a; betaVecContinuous(ii) = pdContinuous.b;
    alphaVecHybrid1(ii) = pdHybrid1.a; betaVecHybrid1(ii) = pdHybrid1.b;
    alphaVecHybrid2(ii) = pdHybrid2.a; betaVecHybrid2(ii) = pdHybrid2.b;
    alphaVecDiscrete(ii) = pdDiscrete.a; betaVecDiscrete(ii) = pdDiscrete.b;
end

% NOTE: Print out the vectors alphaVecContinuous, betaVecContinuous,
%                             alphaVecHybrid1,    betaVecHybrid1,
%                             alphaVecHybrid2,    betaVecHybrid2,
%                             alphaVecDiscrete,   betaVecDiscrete
% in order to fill in the values for cimpval function

fontSize = 15;

width=30;
height=width/5;
figure('paperpositionmode', 'auto', 'units', 'centimeters', 'position', [0 0 width height])

% do the Q-Q plot
pdContinuous = pdObjsContinuous{1};
h1 = subplot(1,3,1); qqplot(rsdmNullDistributionResultsContinuous(:,1), pdContinuous); grid on;
xlabel({['\makebox[4in][c]{Quantiles of ' sprintf('$\\beta(%0.02f, %0.02f)$', ...
         alphaVecContinuous(1), betaVecContinuous(1)) '}'], '\makebox[4in][c]{(a)}'}, ...
    'FontSize', 20, 'Interpreter', 'Latex');
% xlabel(['\makebox[4in][c]{Quantiles of ' sprintf('$\\beta(%0.02f, %0.02f)$', ...
%          alphaVecContinuous(1), betaVecContinuous(1)) '}'], ...
%     'FontSize', 20, 'Interpreter', 'Latex');
ylabel({'Input Samples';'Quantiles'}, 'FontSize', fontSize);
title('M = 100', 'FontSize', fontSize);
h1.FontSize = fontSize;

% pdContinuous = pdObjsContinuous{10};
% h2 = subplot(2,2,2); qqplot(rsdmNullDistributionResultsContinuous(:,10), pdContinuous); grid on;
% xlabel({['\makebox[4in][c]{Quantiles of ' sprintf('$\\beta(%0.02f, %0.02f)$', ...
%     alphaVecContinuous(10), betaVecContinuous(10)) '}'], '\makebox[4in][c]{(b)}'}, ...
%     'FontSize', 20, 'Interpreter', 'Latex');
% ylabel('Quantiles of Input Samples', 'FontSize', fontSize);
% title('M = 1000', 'FontSize', fontSize);
% h2.FontSize = fontSize;

lineWidth = 2; markerSize = 7;

h3 = subplot(1,3,2); 
p3 = plot(M_vec_toPlot, alphaVecContinuous(1:initLen), M_vec_toPlot, alphaVecHybrid1(1:initLen), ...
          M_vec_toPlot, alphaVecHybrid2(1:initLen), M_vec_toPlot, alphaVecDiscrete(1:initLen)); 
h3.XLim = [0,2500];
h3.XTick = [500 1500 2500];
grid on; 
xlabel({'M', '(b)'}, 'FontSize', fontSize); 
% xlabel('M', 'FontSize', fontSize); 
ylabel('\alpha', 'FontSize', fontSize);
% title('(c)', 'FontSize', fontSize);
h3.FontSize = fontSize;
p3(1).LineWidth = lineWidth; p3(1).Marker = 'd'; p3(1).MarkerSize = markerSize;
p3(2).LineWidth = lineWidth; p3(2).Marker = 'v'; p3(2).MarkerSize = markerSize;
p3(3).LineWidth = lineWidth; p3(3).Marker = '*'; p3(3).MarkerSize = markerSize;
p3(4).LineWidth = lineWidth; p3(4).Marker = 'x'; p3(4).MarkerSize = markerSize;

h4 = subplot(1,3,3); 
p4 = plot(M_vec_toPlot, betaVecContinuous(1:initLen), M_vec_toPlot, betaVecHybrid1(1:initLen), ...
          M_vec_toPlot, betaVecHybrid2(1:initLen), M_vec_toPlot, betaVecDiscrete(1:initLen));
h4.XLim = [0,2500];
h4.XTick = [500 1500 2500];
grid on; 
xlabel({'M', '(c)'}, 'FontSize', fontSize); 
% xlabel('M', 'FontSize', fontSize); 
ylabel('\beta', 'FontSize', fontSize);
% title('(d)', 'FontSize', fontSize);
h4.FontSize = fontSize;
p4(1).LineWidth = lineWidth; p4(1).Marker = 'd'; p4(1).MarkerSize = markerSize;
p4(2).LineWidth = lineWidth; p4(2).Marker = 'v'; p4(2).MarkerSize = markerSize;
p4(3).LineWidth = lineWidth; p4(3).Marker = '*'; p4(3).MarkerSize = markerSize;
p4(4).LineWidth = lineWidth; p4(4).Marker = 'x'; p4(4).MarkerSize = markerSize;

legend({'Continuous', 'Hybrid-1', 'Hybrid-2', 'Discrete'}, 'location', 'SouthEast');
