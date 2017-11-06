%% Script which characterizes the CIM algorithm in many different ways
clear;
clc;
close all;
dbstop if error;

masterCfgRun = 1;

% master script configuration so that we can just run the cell and it will
% run all the tests we would like to
runPower_M500              = 0;
runPower_All               = 0;
plotPower_M500_depMeasures = 0;
plotPower_M500_miMeasures  = 0;
plotPower_ss_depMeasures   = 0;
plotPower_ss_miMeasures    = 0;
runPowerSensitivity        = 0;
runPowerAlphaSensitivity   = 0;
plotPowerSensitivity       = 0;
runAlgoSensitivity         = 0;
runAlgoAlphaSensitivity    = 0;
plotAlgoSensitivity        = 0;
runConvergence             = 0;
plotConvergence            = 0;
runPower_test_M500         = 1;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Running Notebook from Master Configuration...\n'),'keepthis','timestamp');

%% Compute the CIM algorithm power vs. other leading measures of dependence and leading measures of DPI
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end

if(~exist('masterCfgRun') || (masterCfgRun==1 && runPower_M500) )
    dispstat(sprintf('Running Power for M=500'),'keepthis','timestamp');
    
    rng(1230);
    
    M = 500;
    minScanIncr = 0.015625;
    mine_c = 15;
    mine_alpha = 0.6;
    rdc_k = 20;
    rdc_s = 1/6;
    knn_1 = 1;
    knn_6 = 6;
    knn_20 = 20;

    nameIdxCorrelationCell = {'CIM', 'dCor','TICe','Corr','RDC','CoS', 'cCor', ...
                              'KNN-1', 'KNN-6', 'KNN-20', 'vME', 'AP'};
    functionHandlesCell = {@cim;
                           @dcor;
                           @mine_interface_tice;
                           @corr;
                           @rdc;
                           @cosdv;
                           @ccor;
                           @KraskovMI_cc_mex;
                           @KraskovMI_cc_mex;
                           @KraskovMI_cc_mex;
                           @vmeMI_interface;
                           @apMI_interface};
    functionArgsCell    = {{minScanIncr};
                           {};
                           {mine_alpha,mine_c,'mic_e'};
                           {};
                           {rdc_k, rdc_s};
                           {};
                           {};
                           {knn_1};
                           {knn_6};
                           {knn_20};
                           {};
                           {};};
    powerCurve = compute_power_curves(M,functionHandlesCell, functionArgsCell);

    % save the data
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_M_%d.mat',M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_M_%d.mat',M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/power_M_%d.mat',M));
    end
end

%% Compute the CIM algorithm power vs. other leading measures of dependence and leading measures of DPI
% for all sample sizes, so we can understand measure's sample size
% requirements
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && runPower_All) )
    dispstat(sprintf('Running Power for M=[25:25:2000]'),'keepthis','timestamp');
    rng(1230);

    MVec = [25:25:2000];
    minScanIncr = 0.015625;
    mine_c = 15;
    mine_alpha = 0.6;
    rdc_k = 20;
    rdc_s = 1/6;
    knn_1 = 1;
    knn_6 = 6;
    knn_20 = 20;

    nameIdxCorrelationCell = {'CIM', 'dCor','TICe','Corr','RDC','CoS', 'cCor', ...
                              'KNN-1', 'KNN-6', 'KNN-20', 'vME', 'AP'};

    functionHandlesCell = {@cim;
                           @dcor;
                           @mine_interface_tice;
                           @corr;
                           @rdc;
                           @cosdv;
                           @ccor;
                           @KraskovMI_cc_mex;
                           @KraskovMI_cc_mex;
                           @KraskovMI_cc_mex;
                           @vmeMI_interface;
                           @apMI_interface};
    functionArgsCell    = {{minScanIncr};
                           {};
                           {mine_alpha,mine_c,'mic_e'};
                           {};
                           {rdc_k, rdc_s};
                           {};
                           {};
                           {knn_1};
                           {knn_6};
                           {knn_20};
                           {};
                           {};};

    num_noise_test_min = 0;
    num_noise_test_max = 30;
    noiseVec = num_noise_test_min:num_noise_test_max;
    numDepMeasures = length(functionHandlesCell);
    numDepTests = 8;
    powerTensor = zeros(length(MVec),numDepMeasures,numDepTests,length(noiseVec));
    for MIdx=1:length(MVec)
        M = MVec(MIdx);
        powerCurve = compute_power_curves(M,functionHandlesCell, functionArgsCell);
        powerTensor(MIdx,:,:,:) = powerCurve;
        % save the data
        if(ispc)
            save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_all.mat');
        elseif(ismac)
            save('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
        else
            save('/home/kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
        end
    end
end

%% Plot the CIM Algorithm Power vs. other leading dependence measures
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && plotPower_M500_depMeasures) )
    % load the data
    if(ispc)
        load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_all.mat');
    elseif(ismac)
        load('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    else
        load('/home/kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    end
    M = 500;
    MIdx = find(MVec==M);

    labels = {'CIM', 'CoS', 'RDC', 'TICe', 'dCor', 'cCor'};
    cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

    num_noise_test_min = 0;
    num_noise_test_max = 20;
    noiseVec = num_noise_test_min:num_noise_test_max;
    powerMat = zeros(length(labels),8,length(noiseVec));

    for labelIdx=1:length(labels)
        label = labels{labelIdx};
        % find which index this corresponds to
        fIdx = find(cellfun(cellfind(label),nameIdxCorrelationCell));
        powerMat(labelIdx,:,:) = squeeze(powerTensor(MIdx,fIdx,:,1:length(noiseVec)));
    end

    noiseVecToPlot = noiseVec/10;

    plotStyle = 1;
    plotPower(powerMat, M, labels, noiseVecToPlot, plotStyle)
end
%% Plot the CIM Algorithm Power vs. other DPI satisfying measures (i.e. measures of Mutual Information)
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && plotPower_M500_miMeasures) )
    % load the data
    if(ispc)
        load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_all.mat');
    elseif(ismac)
        load('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    else
        load('/home/kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    end
    M = 100;
    MIdx = find(MVec==M);

    labels = {'CIM', 'KNN-1', 'KNN-6', 'KNN-20', 'AP', 'vME'};

    num_noise_test_min = 0;
    num_noise_test_max = 20;
    noiseVec = num_noise_test_min:num_noise_test_max;
    powerMat = zeros(length(labels),8,length(noiseVec));

    for labelIdx=1:length(labels)
        label = labels{labelIdx};
        % find which index this corresponds to
        fIdx = find(contains(nameIdxCorrelationCell,label));
        powerMat(labelIdx,:,:) = squeeze(powerTensor(MIdx,fIdx,:,1:length(noiseVec)));
    end

    noiseVecToPlot = noiseVec/10;

    plotStyle = 1;
    plotPower(powerMat, M, labels, noiseVecToPlot, plotStyle)
end
%% Plot the small-sample results for CIM vs. other leading measures of dependence
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && plotPower_ss_depMeasures) )
    % load the data
    if(ispc)
        load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_all.mat');
    elseif(ismac)
        load('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    else
        load('/home/kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    end
    
    labels = {'CIM', 'CoS', 'RDC', 'TICe', 'dCor', 'cCor'};
    cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

    num_noise_test_min = 0;
    num_noise_test_max = 20;
    powerThreshold = 0.7;
    noiseVec = num_noise_test_min:num_noise_test_max;
    sampleSizeAnalysisVec = zeros(length(labels),8,length(noiseVec));
    
    for labelIdx=1:length(labels)
        label = labels{labelIdx};
        % find which index this corresponds to
        fIdx = find(cellfun(cellfind(label),nameIdxCorrelationCell));
        
        for typ=1:numDepTests
            for l=1:length(noiseVec)
                for m=1:length(MVec)
                    M = MVec(m);
                    if(powerTensor(m,fIdx,typ,l)>powerThreshold)
                        % TODO: we can do some interpolation here, so that we
                        % are not restricted to the boundaries of which tests
                        % were run ...
                        sampleSizeAnalysisVec(labelIdx,typ,l) = M;
                        break;
                    end
                end
            end
        end
        
    end

    noiseVecToPlot = noiseVec/10;

    plotStyle = 1;
    plotPower_ss(sampleSizeAnalysisVec, labels, noiseVecToPlot, plotStyle)
end

%% test the power sensitivity of the CIM algorithm to the min-scan-incr
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && runPowerSensitivity) )
    dispstat(sprintf('Testing Power Sensitivity ...'),'keepthis','timestamp');
    rng(1234);
    cimfunc = @cim;
    scanincrsToTest = [0.25, 0.125, .0625, .03125, .015625];
    MVecToTest = 100:100:1000;
    for M=MVecToTest
        powerCurve = cim_power_sensitivity(cimfunc,M,scanincrsToTest);
        % save the results
        if(ispc)
            save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_powerSensitivity_M_%d.mat', fnameStr, M));
        elseif(ismac)
            save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
        else
            save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
        end
    end
end

%% test the power sensitivity of the CIM algorithm to alpha
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && runPowerAlphaSensitivity) )
    dispstat(sprintf('Testing Power Sensitivity ...'),'keepthis','timestamp');
    rng(1234);
    cimfunc = @cim_v2_cc_mex;
    alphasToTest = [0.05 0.1 0.15 0.2 0.25 0.3];
    msiValue = 0.015625;
    MVecToTest = 100:100:1000;
    for M=MVecToTest
        powerCurve = cim_powerAlpha_sensitivity(cimfunc,M,msiValue,alphasToTest);
        % save the results
        if(ispc)
            save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_powerAlphaSensitivity_M_%d.mat', fnameStr, M));
        elseif(ismac)
            save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_powerAlphaSensitivity_M_%d.mat', fnameStr, M));
        else
            save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_powerAlphaSensitivity_M_%d.mat', fnameStr, M));
        end
    end
end


%% Plot the power-sensitivity
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && plotPowerSensitivity) )
    cimVersion = 8;
    MVecToPlot = 100:100:1000;
    num_noise_test_min = 0;
    num_noise_test_max = 20;
    plotPowerSensitivity_withinM(cimVersion,MVecToPlot,num_noise_test_min,num_noise_test_max);
    figtitle(sprintf('Algorithm Power Sensitivity (M=%d - %d)',min(MVecToPlot),max(MVecToPlot)),'FontSize',20);
end


%% Run the algorithm sensitivity analysis
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && runAlgoSensitivity) )
    dispstat(sprintf('Testing Algorithm Sensitivity ...'),'keepthis','timestamp');
    rng(1234);
    scanincrsToTest = [0.25, 0.125, .0625, .03125, .015625];
    cimfunc = @cim;
    MVecToTest = 100:100:1000;
    for M=MVecToTest
        algoSensitivityData = cim_algo_sensitivity(cimfunc,M,scanincrsToTest);    
        % save the results
        if(ispc)
            save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_algoSensitivity_M_%d.mat', fnameStr, M));
        elseif(ismac)
            save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
        else
            save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
        end
    end
end

%% Run the algorithm sensitivity analysis for cim-v2 parametrizing alpha
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && runAlgoAlphaSensitivity) )
    dispstat(sprintf('Testing Algorithm Sensitivity ...'),'keepthis','timestamp');
    rng(1234);
    alphasToTest = [0.05 0.1 0.15 0.2 0.25 0.3];
    msiValue = 0.015625;
    cimfunc = @cim_v2_cc_mex;
    MVecToTest = 100:100:1000;
    for M=MVecToTest
        algoSensitivityData = cim_algoAlpha_sensitivity(cimfunc,M,msiValue,alphasToTest);    
        % save the results
        if(ispc)
            save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_algoSensitivityAlpha_M_%d.mat', fnameStr, M));
        elseif(ismac)
            save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivityAlpha_M_%d.mat', fnameStr, M));
        else
            save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivityAlpha_M_%d.mat', fnameStr, M));
        end
    end
end


%% Plot the Algorithm sensitivity analysis including M
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && plotAlgoSensitivity) )
    cimVersion = 8;
    MVecToPlot = 100:100:1000;
    % plotAlgoSensitivity_acrossM(cimVersion,MVecToPlot);
    num_noise_test_min = 0;
    num_noise_test_max = 20;

    plotAlgoSensitivity_withinM(cimVersion,MVecToPlot,num_noise_test_min,num_noise_test_max);
%     figtitle(sprintf('Algorithm Sensitivity (M=%d - %d)',min(MVecToPlot),max(MVecToPlot)),'FontSize',20);
end
%% simulate the difference between computing monotonic, detecting the region, and theoretical CIM value
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && runConvergence) )
    dispstat(sprintf('Testing Convergence ...'),'keepthis','timestamp');
    rng(12345);
    cimfunc = @cim;
    MVec = [100:100:1000 2000 5000 10000];
    fArgs = {0.015625};  % the minScanIncr value

    for M=MVec
        [linearDep,quadraticDep,cubicDep,...
              sinusoidalDep,hiFreqSinDep,fourthRootDep,...
              circleDep,stepDep,indep] = ...
            compute_cim_convergence(cimfunc, fArgs, M);
        % save the data
        if(ispc)
            save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_convergence_M_%d.mat', M));
        elseif(ismac)
            save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_convergence_M_%d.mat', M));
        else
            save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_convergence_M_%d.mat', M));
        end
    end
end
%% Plot the results
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end
if(~exist('masterCfgRun') || (masterCfgRun==1 && plotConvergence) )
    
    plotValues = 0;  % if 1, we plot the values, if 0, we plot the bias
    
    MVecDataAvailable = [100:100:1000 2000 5000];  % add 5000 and 10000 to this as they come online
    tol = 0.015;
    noiseMinToAnalyze = 0; noiseMaxToAnalyze = 20;
    noiseVecToAnalyze = noiseMinToAnalyze:noiseMaxToAnalyze;  

    % do a +1 to account for matlab indexing
    noiseVecToAnalyze = noiseVecToAnalyze + 1;
    aggregatefunc = @mean;

    numDep = 8;
    depTestVec = ones(1,numDep);
    MVecResults = zeros(1,numDep);
    for MIdx=1:length(MVecDataAvailable)
        M = MVecDataAvailable(MIdx);
        if(ispc)
            load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_convergence_M_%d.mat', M));
        elseif(ismac)
            load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_convergence_M_%d.mat', M));
        else
            load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_convergence_M_%d.mat', M));
        end

        % compute error between theoretical CIM and actual CIM for each of the
        % dependencies, if the total error is within the sseTol value, then we
        % store that as the M Value for which the error is acceptable.

        if(depTestVec(1))
    %         linearDepErr = (linearDep(2,noiseToTest)-linearDep(3,noiseToTest)).^2;
            linearDepErr = abs(linearDep(2,noiseVecToAnalyze)-linearDep(3,noiseVecToAnalyze));
            if(aggregatefunc(linearDepErr)<=tol)
                MVecResults(1) = M;
                depTestVec(1) = 0;
                linearDepToPlot = linearDep;
            end
        end

        if(depTestVec(2))
    %         quadraticDepErr = (quadraticDep(2,noiseToTest)-quadraticDep(3,noiseToTest)).^2;
            quadraticDepErr = abs(quadraticDep(2,noiseVecToAnalyze)-quadraticDep(3,noiseVecToAnalyze));
            if(aggregatefunc(quadraticDepErr)<=tol)
                MVecResults(2) = M;
                depTestVec(2) = 0;
                quadraticDepToPlot = quadraticDep;
            end
        end

        if(depTestVec(3))
    %         cubicDepErr = (cubicDep(2,noiseToTest)-cubicDep(3,noiseToTest)).^2;
            cubicDepErr = abs(cubicDep(2,noiseVecToAnalyze)-cubicDep(3,noiseVecToAnalyze));
            if(aggregatefunc(cubicDepErr)<=tol)
                MVecResults(3) = M;
                depTestVec(3) = 0;
                cubicDepToPlot = cubicDep;
            end
        end

        if(depTestVec(4))
    %         sinusoidalDepErr = (sinusoidalDep(2,noiseToTest)-sinusoidalDep(3,noiseToTest)).^2;
            sinusoidalDepErr = abs(sinusoidalDep(2,noiseVecToAnalyze)-sinusoidalDep(3,noiseVecToAnalyze));
            if(aggregatefunc(sinusoidalDepErr)<=tol)
                MVecResults(4) = M;
                depTestVec(4) = 0;
                sinusoidalDepToPlot = sinusoidalDep;
            end
        end

        if(depTestVec(5))
    %         hiFreqSinDepErr = (hiFreqSinDep(2,noiseToTest)-hiFreqSinDep(3,noiseToTest)).^2;
    %         hiFreqSinDep(2,1) = 1;
            hiFreqSinDepErr = abs(hiFreqSinDep(2,noiseVecToAnalyze)-hiFreqSinDep(3,noiseVecToAnalyze));
            if(aggregatefunc(hiFreqSinDepErr)<=tol)
                MVecResults(5) = M;
                depTestVec(5) = 0;
                hiFreqSinDepToPlot = hiFreqSinDep;
            end
        end

        if(depTestVec(6))
    %         fourthRootDepErr = (fourthRootDep(2,noiseToTest)-fourthRootDep(3,noiseToTest)).^2;
            fourthRootDepErr = abs(fourthRootDep(2,noiseVecToAnalyze)-fourthRootDep(3,noiseVecToAnalyze));
            if(aggregatefunc(fourthRootDepErr)<=tol)
                MVecResults(6) = M;
                depTestVec(6) = 0;
                fourthRootDepToPlot = fourthRootDep;
            end
        end

        if(depTestVec(7))
    %         circleDepErr = (circleDep(2,noiseToTest)-circleDep(3,noiseToTest)).^2;
            circleDepErr = abs(circleDep(2,noiseVecToAnalyze)-circleDep(3,noiseVecToAnalyze));
            if(aggregatefunc(circleDepErr)<=.01)
                MVecResults(7) = M;
                depTestVec(7) = 0;
                circleDepToPlot = circleDep;
            end
        end

        if(depTestVec(8))
    %         stepDepErr = (stepDep(2,noiseToTest)-stepDep(3,noiseToTest)).^2;
            stepDep(2,1) = 1; % due to a bug in the way we simulated we need to force
                              % this value to 1
            stepDepErr = abs(stepDep(2,noiseVecToAnalyze)-stepDep(3,noiseVecToAnalyze));
            if(aggregatefunc(stepDepErr)<=tol)
                MVecResults(8) = M;
                depTestVec(8) = 0;
                stepDepToPlot = stepDep;
            end
        end
    end

    % use max M for any remainders that didn't meet the SSE
    for depTestValIdx=1:length(depTestVec)
        if(depTestVec(depTestValIdx)==1)
            M = max(MVecDataAvailable);
%             if(ispc)
%                 load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_fourier_convergence_%s_M_%d.mat', fnameStr, M));
%             elseif(ismac)
%                 load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_fourier_convergence_%s_M_%d.mat', fnameStr, M));
%             else
%                 load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_fourier_convergence_%s_M_%d.mat', fnameStr, M));
%             end
%             if(ispc)
%                 load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_convergence_%s_M_%d.mat', fnameStr, M));
%             elseif(ismac)
%                 load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_convergence_%s_M_%d.mat', fnameStr, M));
%             else
%                 load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_convergence_%s_M_%d.mat', fnameStr, M));
%             end
            if(ispc)
                load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_convergence_M_%d.mat', M));
            elseif(ismac)
                load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_convergence_M_%d.mat', M));
            else
                load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_convergence_M_%d.mat', M));
            end
            
            MVecResults(depTestValIdx) = M;
            switch depTestValIdx
                case 1
                    linearDepToPlot = linearDep;
                case 2
                    quadraticDepToPlot = quadraticDep;
                case 3
                    cubicDepToPlot = cubicDep;
                case 4
                    sinusoidalDepToPlot = sinusoidalDep;
                case 5
                    hiFreqSinDepToPlot = hiFreqSinDep;
                case 6
                    fourthRootDepToPlot = fourthRootDep;
                case 7
                    circleDepToPlot = circleDep;
                case 8
                    stepDepToPlot = stepDep;
                otherwise
                    error('UNK!');
            end
        end
    end

    M_inlet = 200;
    inletX = linspace(0,1,M_inlet);
    inletT = linspace(0,2*pi,M_inlet);
    inletData = zeros(8,M_inlet);
    inletData(1,:) = inletX;
    inletData(2,:) = 4*(inletX-.5).^2;
    inletData(3,:) = 128*(inletX-1/3).^3-48*(inletX-1/3).^3-12*(inletX-1/3);
    inletData(4,:) = sin(4*pi*inletX);
    inletData(5,:) = sin(16*pi*inletX);
    inletData(6,:) = inletX.^(1/4);
    inletData(8,:) = (inletX > 0.5);
    inset_bufX = .097; inset_bufY = 0.28;
    inset_width = 0.06; inset_height = 0.06;

    noiseVecPlot = noiseVecToAnalyze-1;
    lineWidthVal = 3;

    figure;
    h = subplot(2,4,1);
    % hh1 = plot(noiseToTest,linearDepToPlot(1,noiseToTest),'o-.', ...
    %      noiseToTest,linearDepToPlot(2,noiseToTest),'+-.', ...
    %      noiseToTest,linearDepToPlot(3,noiseToTest),'d-.');
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,linearDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,linearDepToPlot(3,noiseVecToAnalyze),'d-.');
        hLegend = legend('CIM','$$\widehat{CIM}$$');
        set(hLegend,'Interpreter','latex','Location','SouthWest')
    else
        hh1 = plot(noiseVecPlot/10,abs(linearDepToPlot(2,noiseVecToAnalyze)-...
                                       linearDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(1)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    h.FontSize = 20;
    inletIdx = 1;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
    ax.XLim = [min(inletX) max(inletX)];
    ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

    h = subplot(2,4,2);
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,quadraticDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,quadraticDepToPlot(3,noiseVecToAnalyze),'d-.');
    else
        hh1 = plot(noiseVecPlot/10,abs(quadraticDepToPlot(2,noiseVecToAnalyze)-...
                                       quadraticDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(2)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    % hh1(3).LineWidth = 1.5; 
    h.FontSize = 20;
    inletIdx = 2;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
    ax.XLim = [min(inletX) max(inletX)];
    ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

    h = subplot(2,4,3);
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,cubicDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,cubicDepToPlot(3,noiseVecToAnalyze),'d-.');
    else
        hh1 = plot(noiseVecPlot/10,abs(cubicDepToPlot(2,noiseVecToAnalyze)-...
                                       cubicDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(3)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    % hh1(3).LineWidth = 1.5; 
    h.FontSize = 20;
    inletIdx = 3;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
    ax.XLim = [min(inletX) max(inletX)];
    ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

    h = subplot(2,4,4);
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,sinusoidalDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,sinusoidalDepToPlot(3,noiseVecToAnalyze),'d-.');
    else
        hh1 = plot(noiseVecPlot/10,abs(sinusoidalDepToPlot(2,noiseVecToAnalyze)-...
                                       sinusoidalDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(4)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    % hh1(3).LineWidth = 1.5; 
    h.FontSize = 20;
    inletIdx = 4;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
    ax.XLim = [min(inletX) max(inletX)];
    ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

    h = subplot(2,4,5);
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,hiFreqSinDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,hiFreqSinDepToPlot(3,noiseVecToAnalyze),'d-.');
    else
        hh1 = plot(noiseVecPlot/10,abs(hiFreqSinDepToPlot(2,noiseVecToAnalyze)-...
                                       hiFreqSinDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(5)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    % hh1(3).LineWidth = 1.5; 
    h.FontSize = 20;
    inletIdx = 5;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
    ax.XLim = [min(inletX) max(inletX)];
    ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

    h = subplot(2,4,6);
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,fourthRootDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,fourthRootDepToPlot(3,noiseVecToAnalyze),'d-.');
    else
        hh1 = plot(noiseVecPlot/10,abs(fourthRootDepToPlot(2,noiseVecToAnalyze)-...
                                       fourthRootDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(6)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    % hh1(3).LineWidth = 1.5; 
    h.FontSize = 20;
    inletIdx = 6;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
    ax.XLim = [min(inletX) max(inletX)];
    ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

    h = subplot(2,4,7);
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,circleDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,circleDepToPlot(3,noiseVecToAnalyze),'d-.');
    else
        hh1 = plot(noiseVecPlot/10,abs(circleDepToPlot(2,noiseVecToAnalyze)-...
                                       circleDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(7)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    % hh1(3).LineWidth = 1.5; 
    h.FontSize = 20;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(cos(inletT),sin(inletT), 'k', 'LineWidth', 2);
    ax.XLim = [min(cos(inletT)) max(cos(inletT))];
    ax.YLim = [min(sin(inletT)) max(sin(inletT))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

    h = subplot(2,4,8);
    if(plotValues)
        hh1 = plot(noiseVecPlot/10,stepDepToPlot(2,noiseVecToAnalyze),'+-.', ...
                   noiseVecPlot/10,stepDepToPlot(3,noiseVecToAnalyze),'d-.');
    else
        hh1 = plot(noiseVecPlot/10,abs(stepDepToPlot(2,noiseVecToAnalyze)-...
                                       stepDepToPlot(3,noiseVecToAnalyze)),'+-.');
    end
    grid on;
%     xlabel('Noise','FontSize',20);
    title(sprintf('min(M)=%d',MVecResults(8)),'FontSize',20);
    hh1(1).LineWidth = lineWidthVal; 
    if(plotValues)
        hh1(2).LineWidth = lineWidthVal; 
    end
    % hh1(3).LineWidth = 1.5; 
    h.FontSize = 20;
    inletIdx = 8;
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
    ax.XLim = [min(inletX) max(inletX)];
    ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];
    
    [~,hL] = suplabel('Noise','x');
    set(hL,'FontSize',20);
    if(plotValues)
        [~,hL] = suplabel('','y');
        hL.FontSize = 20;
        hL.Interpreter = 'Latex';
    else
        [~,hL] = suplabel('Bias','y');
        hL.FontSize = 20;
    end

    % subplot(3,3,9);
    % hh1 = plot(noiseVec,indep(1,:),'o-.', ...
    %      noiseVec,indep(2,:),'+-.', ...
    %      noiseVec,indep(3,:),'d-.');
    % grid on;
    % xlabel('noise');
    % title('Independence');
    % hh1(1).LineWidth = 1.5; 
    % hh1(2).LineWidth = 1.5; 
    % hh1(3).LineWidth = 1.5; 
end

%% Test CIM heuristic vs CIM 
if(~exist('masterCfgRun'))  % means we are running the cell independently
    clear;
    clc;
    close all;
    dbstop if error;
    dispstat('','init'); % One time only initialization
end

if(~exist('masterCfgRun') || (masterCfgRun==1 && runPower_test_M500) )
    dispstat(sprintf('Running Power for M=500'),'keepthis','timestamp');
    
    rng(1230);
    
    M = 500;
    alpha = 0.2;
    minScanIncr = 0.015625;
    nsim_null = 300;
    nsim_alt = 300;
    num_noise_test_min = 0;
    num_noise_test_max = 20;
    
    nameIdxCorrelationCell = {'CIM', 'CIMv2'};
    functionHandlesCell = {@cim;
                           @cim_v2_cc_mex;};
    functionArgsCell    = {{minScanIncr};
                           {minScanIncr, alpha};
                          };
    powerCurve = compute_power_curves(M,functionHandlesCell, functionArgsCell,...
                                      nsim_null,nsim_alt,...
                                      num_noise_test_min,num_noise_test_max);

    % save the data
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_comparison_power_alpha0.3_M_%d.mat',M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_comparison_power_alpha0.3_M_%d.mat',M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_comparison_power_alpha0.3_M_%d.mat',M));
    end
end

%% Plot the CIM Algorithm Power vs. other leading dependence measures
% load the data

clear;
clc;
M = 500;
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_comparison_power_alpha0.3_M_%d.mat',M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_comparison_power_alpha0.3_M_%d.mat',M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_comparison_power_alpha0.3_M_%d.mat',M));
end

labels = {'CIM', 'CIMv2'};
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

num_noise_test_min = 0;
num_noise_test_max = 20;
noiseVec = num_noise_test_min:num_noise_test_max;

noiseVecToPlot = noiseVec/10;

plotStyle = 1;
plotPower(powerCurve, M, labels, noiseVecToPlot, plotStyle)