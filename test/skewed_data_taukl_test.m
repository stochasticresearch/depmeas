%% Test how tau-kl is affected by skewed distributions for continuous vs. binary
% Note, we only care about binary here, so we generate bernoulli
% distribution, but in the future, we also want to simulate non-binary
% skewed output variables

clear;
clc;
dbstop if error;

scenarios = {'left-skew','no-skew','right-skew'};
tauVec = linspace(0.01,0.99,20);
copulas = {'Gaussian','Frank','Gumbel','Clayton'};
M = 500;
numMCSims = 100;

% manually generate a left-skewed and right-skewed data, from which we
% construct an empirical cdf
leftSkewData = pearsrnd(0,1,-1,3,5000,1);
rightSkewData = pearsrnd(0,1,1,3,5000,1);

[fLeftSkew,xiLeftSkew] = emppdf(leftSkewData,0);
FLeftSkew = empcdf(xiLeftSkew,0);
leftSkewContinuousDistInfo = rvEmpiricalInfo(xiLeftSkew,fLeftSkew,FLeftSkew,0);
[fRightSkew,xiRightSkew] = emppdf(rightSkewData,0);
FRightSkew = empcdf(xiRightSkew,0);
rightSkewContinuousDistInfo = rvEmpiricalInfo(xiRightSkew,fRightSkew,FRightSkew,0);

resVecTauB  = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecTauN  = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecTauKL = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecCIM   = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecKNN1  = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecKNN6  = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecKNN20 = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecVME   = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecAP    = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecEntropyMI= zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));

msi = 0.015625; alpha = 0.2;

sampleDataMat = zeros(M,2,length(copulas),length(tauVec),length(scenarios),length(scenarios));

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

aa = 1; 
for continuousDistScenario=scenarios
    bb = 1; 
    for discreteDistScenario=scenarios
        cc = 1; 
        for tau=tauVec
            dd = 1; 
            for cop=copulas
                dispstat(sprintf('Simulating for {%s,%s} tau=%0.02f >> %s',cell2str(continuousDistScenario),cell2str(discreteDistScenario),tau,cell2str(cop)), ...
                         'keepthis', 'timestamp');
        
                iTau = copulaparam(cop,tau);
                parfor mcSimNum=1:numMCSims
%                 for mcSimNum=1:numMCSims
                    % generate U
                    U = copularnd(cop,iTau,M);
                    
                    % generate F_X
                    if(strcmpi('left-skew',continuousDistScenario))
                        X = zeros(length(U),1);
                        % TODO: vectorize
                        for ii=1:length(U)
                            X(ii) = leftSkewContinuousDistInfo.icdf(U(ii,1));
                        end
                    elseif(strcmpi('no-skew',continuousDistScenario))
                        distObj = makedist('Normal');
                        X = icdf(distObj,U(:,1));
                    elseif(strcmpi('right-skew',continuousDistScenario))
                        X = zeros(length(U),1);
                        % TODO: vectorize
                        for ii=1:length(U)
                            X(ii) = rightSkewContinuousDistInfo.icdf(U(ii,1));
                        end
                    end
                    
                    % generate F_Y
                    if(strcmpi('left-skew',discreteDistScenario))
                        distObj = makedist('Multinomial','probabilities',[0.1,0.9]);
                    elseif(strcmpi('no-skew',discreteDistScenario))
                        distObj = makedist('Multinomial','probabilities',[0.5,0.5]);
                    elseif(strcmpi('right-skew',discreteDistScenario))
                        distObj = makedist('Multinomial','probabilities',[0.9,0.1]);
                    end
                    Y = icdf(distObj,U(:,2));
                    [X,Y] = pobs_sorted_cc(X,Y);

                    % compute tau, tau_kl, tau_N and record
                    resVecTauN(mcSimNum,dd,cc,bb,aa) = corr(X,Y,'type','kendall');
                    resVecTauB(mcSimNum,dd,cc,bb,aa) = ktaub([X Y], 0.05, 0);
                    resVecTauKL(mcSimNum,dd,cc,bb,aa) = taukl_cc(X,Y,0,1,0);
                    resVecCIM(mcSimNum,dd,cc,bb,aa) = cim_cc_mex(X,Y,msi,alpha,0,1,0);
                    
                    resVecKNN1(mcSimNum,dd,cc,bb,aa) = KraskovMI_cc_mex(X,Y,1);
                    resVecKNN6(mcSimNum,dd,cc,bb,aa) = KraskovMI_cc_mex(X,Y,6);
                    resVecKNN20(mcSimNum,dd,cc,bb,aa) = KraskovMI_cc_mex(X,Y,20);
                    resVecVME(mcSimNum,dd,cc,bb,aa) = vmeMI_interface(X,Y);
                    resVecAP(mcSimNum,dd,cc,bb,aa) = apMI_interface(X,Y);

                    X_continuous = 1;
                    resVecEntropyMI(mcSimNum,dd,cc,bb,aa) = h_mi_interface(X,Y,X_continuous);
                end                
                
                %%%%%%%%%%%%%%%%%%%% MESSY CODE !!!!!! %%%%%%%%%%%%%%%%%%%%
                % COPIED FROM ABOVE, b/c parfor suks sometimes :/
                % generate U
                U = copularnd(cop,iTau,M);

                % generate F_X
                if(strcmpi('left-skew',continuousDistScenario))
                    X = zeros(length(U),1);
                    % TODO: vectorize
                    for ii=1:length(U)
                        X(ii) = leftSkewContinuousDistInfo.icdf(U(ii,1));
                    end
                elseif(strcmpi('no-skew',continuousDistScenario))
                    distObj = makedist('Normal');
                    X = icdf(distObj,U(:,1));
                elseif(strcmpi('right-skew',continuousDistScenario))
                    X = zeros(length(U),1);
                    % TODO: vectorize
                    for ii=1:length(U)
                        X(ii) = rightSkewContinuousDistInfo.icdf(U(ii,1));
                    end
                end

                % generate F_Y
                numIndepTrials = 1;
                if(strcmpi('left-skew',discreteDistScenario))
                    distObj = makedist('Multinomial','probabilities',[0.1,0.9]);
                elseif(strcmpi('no-skew',discreteDistScenario))
                    distObj = makedist('Multinomial','probabilities',[0.5,0.5]);
                elseif(strcmpi('right-skew',discreteDistScenario))
                    distObj = makedist('Multinomial','probabilities',[0.9,0.1]);
                end
                Y = icdf(distObj,U(:,2));
                sampleDataMat(:,1,dd,cc,bb,aa) = X;
                sampleDataMat(:,2,dd,cc,bb,aa) = Y;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dd = dd + 1;
            end
            cc = cc + 1;
        end
        bb = bb + 1;
    end
    aa = aa + 1;
end

% save the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\skewed_data\\binary_output_class_entropybased.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class_entropybased.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class_entropybased.mat');
end

%% plot the bias for tau-variants
clear;
clc;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\skewed_data\\binary_output_class.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
end


sampleData1_subplotIdxVec = [1,3,5,   13,15,17, 25,27,29];
sampleData2_subplotIdxVec = [2,4,6,   14,16,18, 26,28,30];
sampleData3_subplotIdxVec = [7,9,11,  19,21,23, 31,33,35];
results_subplotIdxVec     = [8,10,12, 20,22,24, 32,34,36];

for ii=1:length(copulas)
    copToVis = copulas{ii};
    % find the appropriate indices
    dd = find(contains(copulas,copToVis));

    f = figure;
    subplotIdx = 1;
    for aa=1:length(scenarios)
        for bb=1:length(scenarios)
            tauNVec = mean(squeeze(resVecTauN(:,dd,:,bb,aa)));
            tauBVec = mean(squeeze(resVecTauB(:,dd,:,bb,aa)));
            tauKLVec = mean(squeeze(resVecTauKL(:,dd,:,bb,aa)));
            cimVec = mean(squeeze(resVecCIM(:,dd,:,bb,aa)));
            
            % plot the first sample data
            hh = subplot(6,6,sampleData1_subplotIdxVec(subplotIdx));
            tauIdx = 1;
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
            
            % plot the second sample data
            hh = subplot(6,6,sampleData2_subplotIdxVec(subplotIdx));
            tauIdx = floor(length(tauVec)/2);
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
            
            % plot the third sample data
            hh = subplot(6,6,sampleData3_subplotIdxVec(subplotIdx));
            tauIdx = length(tauVec);
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
            
            % plot the results
            subplot(6,6,results_subplotIdxVec(subplotIdx));
            plot(tauVec,tauBVec,tauVec,tauNVec,tauVec,tauKLVec,tauVec,cimVec,'+');
            xlabel('\tau');
            title(sprintf('%s | %s', scenarios{bb}, scenarios{aa}));
            if(subplotIdx==9)
%                 legend('\tau_b','\tau_N','\tau_{KL}', 'location', 'northwest');
                legend('\tau_b','\tau_N','\tau_{KL}','CIM', 'location', 'northwest');
            end
            grid on;
            
            subplotIdx = subplotIdx + 1;
        end
    end
    tightfig;
    set(f, 'Name', copToVis);
end

%% plot the bias for MI-variants
clear;
clc;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\skewed_data\\binary_output_class.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
end

sampleData1_subplotIdxVec = [1,3,5,   13,15,17, 25,27,29];
sampleData2_subplotIdxVec = [2,4,6,   14,16,18, 26,28,30];
sampleData3_subplotIdxVec = [7,9,11,  19,21,23, 31,33,35];
results_subplotIdxVec     = [8,10,12, 20,22,24, 32,34,36];

for ii=1:length(copulas)
    copToVis = copulas{ii};
    % find the appropriate indices
    dd = find(contains(copulas,copToVis));

    f = figure;
    subplotIdx = 1;
    for aa=1:length(scenarios)
        for bb=1:length(scenarios)
            resVecKNN1Vec = mean(squeeze(resVecKNN1(:,dd,:,bb,aa)));
            resVecKNN6Vec = mean(squeeze(resVecKNN6(:,dd,:,bb,aa)));
            resVecKNN20Vec = mean(squeeze(resVecKNN20(:,dd,:,bb,aa)));
            resVecVMEVec = mean(squeeze(resVecVME(:,dd,:,bb,aa)));
            resVecAPVec = mean(squeeze(resVecAP(:,dd,:,bb,aa)));
            
            % plot the first sample data
            hh = subplot(6,6,sampleData1_subplotIdxVec(subplotIdx));
            tauIdx = 1;
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
            
            % plot the second sample data
            hh = subplot(6,6,sampleData2_subplotIdxVec(subplotIdx));
            tauIdx = floor(length(tauVec)/2);
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
            
            % plot the third sample data
            hh = subplot(6,6,sampleData3_subplotIdxVec(subplotIdx));
            tauIdx = length(tauVec);
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
            
            % plot the results
            subplot(6,6,results_subplotIdxVec(subplotIdx));
            plot(tauVec,resVecKNN1Vec,tauVec,resVecKNN6Vec,tauVec,resVecKNN20Vec, ...
                 tauVec,resVecVMEVec,tauVec,resVecAPVec);
            xlabel('\tau');
            title(sprintf('%s | %s', scenarios{bb}, scenarios{aa}));
            if(subplotIdx==9)
                legend('KNN-1','KNN-6','KNN-20', 'vME', 'AP', 'location', 'northwest');
            end
            grid on;
            
            subplotIdx = subplotIdx + 1;
        end
    end
    tightfig;
    set(f, 'Name', copToVis);
end

%% Figure 1 - SHow SKewed Data nature
clear;
clc;
close all;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\skewed_data\\binary_output_class.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
end

copToVis = 'Gaussian';
% find the appropriate indices
dd = find(contains(copulas,copToVis));
            
f = figure;

subplotIdx=1;
for scenarioIdx=1:length(scenarios)
    skewIdx = find(contains(scenarios,scenarios{scenarioIdx}));
    % plot the first sample data
    hh = subplot(3,3,subplotIdx);
    tauIdx = 1;
    X = sampleDataMat(:,1,dd,tauIdx,skewIdx,skewIdx);
    Y = sampleDataMat(:,2,dd,tauIdx,skewIdx,skewIdx);
    scatter(X,Y);
    set(hh,'XTick',[]); set(hh,'YTick',[]); 
    axis([min(X),max(X),0.8,2.2]);
    if(subplotIdx==1 || subplotIdx==4 || subplotIdx==7)
        ylabel(scenarios{scenarioIdx});
    end
    if(subplotIdx==7)
        xlabel(sprintf('\\tau=%0.02f',tauVec(tauIdx)));
    end
    subplotIdx = subplotIdx + 1;

    % plot the second sample data
    hh = subplot(3,3,subplotIdx);
    tauIdx = floor(length(tauVec)/2);
    X = sampleDataMat(:,1,dd,tauIdx,skewIdx,skewIdx);
    Y = sampleDataMat(:,2,dd,tauIdx,skewIdx,skewIdx);
    scatter(X,Y);
    set(hh,'XTick',[]); set(hh,'YTick',[]); 
    axis([min(X),max(X),0.8,2.2]);
    if(subplotIdx==8)
        xlabel(sprintf('\\tau=%0.02f',tauVec(tauIdx)));
    end
    subplotIdx = subplotIdx + 1;
    
    % plot the third sample data
    hh = subplot(3,3,subplotIdx);
    tauIdx = length(tauVec);
    X = sampleDataMat(:,1,dd,tauIdx,skewIdx,skewIdx);
    Y = sampleDataMat(:,2,dd,tauIdx,skewIdx,skewIdx);
    scatter(X,Y);
    set(hh,'XTick',[]); set(hh,'YTick',[]); 
    axis([min(X),max(X),0.8,2.2]);
    if(subplotIdx==9)
        xlabel(sprintf('\\tau=%0.02f',tauVec(tauIdx)));
    end
    subplotIdx = subplotIdx + 1;
    
end
tightfig;
% suptitle(copToVis);

%% make the plots for the paper - Figure 2 - synth data results
clear;
clc;
close all;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\skewed_data\\binary_output_class.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
end

fontSize = 16;

for dd=1:length(copulas)
% for dd=1:1
    copToVis = copulas{dd};
    width = 25; height = width/3.5;
    figure('paperpositionmode', 'auto', 'units', 'centimeters', 'position', [0 0 width height])
    subplotIdx=1;
    linkList = [];
    for scenarioIdx=1:length(scenarios)
        skewIdx = find(contains(scenarios,scenarios{scenarioIdx}));

        % plot the first sample data
        hh = subplot(1,3,subplotIdx); linkList = [linkList hh];
        tauNVec = mean(squeeze(resVecTauN(:,dd,:,skewIdx,skewIdx)));
        cimVec = mean(squeeze(resVecCIM(:,dd,:,skewIdx,skewIdx))); 
        resVecKNN20Vec = mean(squeeze(resVecKNN20(:,dd,:,skewIdx,skewIdx)));
        resVecVMEVec = mean(squeeze(resVecVME(:,dd,:,skewIdx,skewIdx)));
        resVecAPVec = mean(squeeze(resVecAP(:,dd,:,skewIdx,skewIdx)));
        resEntropyMIVec = mean(squeeze(resVecEntropyMI(:,dd,:,skewIdx,skewIdx)));
        
        hln(:,1) = plot(tauVec,cimVec, tauVec,tauNVec,tauVec,resVecKNN20Vec, ...
                     tauVec,resVecVMEVec,tauVec,resVecAPVec, ...
                     tauVec,resEntropyMIVec, ...
                     tauVec,tauVec,'*k');
        hln(1).LineWidth = 2.5;
        axis([0 1 0 1]);
        
        grid on;
        if(subplotIdx==1)
            ylabel('Measure','FontSize',fontSize);
        end
        % remove y ticks when plotting subplot 2 & 3
        if(subplotIdx==2 || subplotIdx==3)
            hh.YTickLabel = [];
        end
        if(subplotIdx==1)
            title('Left-Skew','FontSize',fontSize);
        elseif(subplotIdx==2)
%             title({copToVis,'No-Skew'},'FontSize',fontSize);
            title('No-Skew','FontSize',fontSize);
        elseif(subplotIdx==3)
            title('Right-Skew','FontSize',fontSize);
        end
        if(subplotIdx==2)
            xlabel('\tau','FontSize',fontSize);            
            if(dd==1)
                legendCell = {'CIM','\tau_N','KNN_{20}','vME','AP','H'};
                [hl(1).leg, hl(1).obj, hl(1).hout, hl(1).mout] = ...
                    legendflex(hln(:,1), legendCell, 'anchor', {'nw','nw'}, ...
                    'buffer', [10 30], ...
                    'ncol', 2, ...
                    'fontsize', fontSize-3, ...
                    'xscale', 0.4, ...
                    'box', 'off');
            end
        end
        
        subplotIdx = subplotIdx + 1;
    end
    linkaxes(linkList,'xy');
%     h = suptitle(copToVis);
%     h.FontSize = fontSize+2;
end

%% Compare the skewed hybrid data versus the entropy & conditional entropy based methods of estimating MI


%% compare taukl to CIM
clear;
clc;
close all;

tau = 0.6;
cop = 'Gaussian';
iTau = copulaparam(cop,tau);
M = 500;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

numMCSim = 1000;
tau_vec = zeros(1,numMCSim);
cim_vec = zeros(1,numMCSim);
probVecCell = {[0.05,0.95],[0.2,0.8],[0.5,0.5],[0.8,0.2],[0.95,0.05]};
for jj=1:length(probVecCell)
    probVec = probVecCell{jj};
    for ii=1:numMCSim
        dispstat(sprintf('%d/%d',ii,numMCSim),'timestamp');
        U = copularnd(cop,iTau,M);
        distObj = makedist('Normal');
        X = icdf(distObj,U(:,1));                  
        distObj = makedist('Multinomial','probabilities',probVec);
        Y = icdf(distObj,U(:,2));
        [u,v] = pobs_sorted_cc(X,Y);
        [v_reverse_sorted,u_reverse_sorted] = pobs_sorted_cc(Y,X);

        msi = 0.015625;
        alpha = 0.2;

        tau_val = taukl_cc(u,v,0,1,0);
        cim_val = cim_cc_mex(u,v,msi,alpha,0,1,0);
        correct = abs(tau_val - cim_val)<=0.01;
        if (~correct)
            dispstat(sprintf('!!! tau=%0.02f cim=%0.02f\n',tau_val,cim_val),'keepthis','timestamp');
            pause;
        end
        tau_vec(ii) = tau_val;
        cim_vec(ii) = cim_val;
    end
    total_err = sum(tau_vec-cim_vec);
    dispstat(sprintf('err{[%0.02f,%0.02f]}=%0.02f',probVec(1),probVec(2),total_err),'keepthis','timestamp');
end

%%
%% compare taukl_cc to taukl
clear;
clc;

tau = 0.6;
cop = 'Gaussian';
iTau = copulaparam(cop,tau);
M = 500;

numMCSim = 5000;

sumCorrect = 0;
for simNum=1:numMCSim
    U = copularnd(cop,iTau,M);
    distObj = makedist('Normal');
    X = icdf(distObj,U(:,1));                  
    distObj = makedist('Multinomial','probabilities',[0.05,0.95]);
    Y = icdf(distObj,U(:,2));

    % intentionally swap u and v
    [u,v] = pobs_sorted_cc(Y,X);

    % we have to tell it Y is the continuous variable, b/c we swapped X and Y in the pobs computation
    % (this is the last argument of 1 instead of 0)
    ccval = taukl_cc(u,v,0,1,1);
    regval = taukl(X,Y);

    correct = abs(ccval - regval)<=eps;
    sumCorrect = sumCorrect + correct;
end
sumCorrect