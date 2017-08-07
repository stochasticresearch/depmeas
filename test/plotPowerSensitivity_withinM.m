function [] = plotPowerSensitivity_withinM(cimVersion,MVecToPlot,num_noise_test_min,num_noise_test_max)

switch cimVersion
    case 1
        fnameStr = 'cim';
    case 3
        fnameStr = 'cimv3';
    case 4
        fnameStr = 'cimv4';
    case 5
        fnameStr = 'cimv5';
    case 6
        fnameStr = 'cimv6';
    case 7
        fnameStr = 'cimv7';
    case 8
        fnameStr = 'cimv8';
    otherwise
end

linearDepCell = cell(1,length(MVecToPlot));
quadraticDepCell = cell(1,length(MVecToPlot));
cubicDepCell = cell(1,length(MVecToPlot));
sinusoidalDepCell = cell(1,length(MVecToPlot));
hiFreqSinDepCell = cell(1,length(MVecToPlot));
fourthRootDepCell = cell(1,length(MVecToPlot));
circleDepCell = cell(1,length(MVecToPlot));
stepDepCell = cell(1,length(MVecToPlot));

numLoopsRun = 1;
for MIdx=1:length(MVecToPlot)
    M = MVecToPlot(MIdx);
    % load the data
    if(ispc)
        load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_powerSensitivity_M_%d.mat', fnameStr, M));
    elseif(ismac)
        load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
    else
        load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
    end
    noiseVec = num_noise_test_min:num_noise_test_max;

    linearDep = zeros(2,length(noiseVec));      % 1 - min
                                                % 2 - max
    quadraticDep = zeros(2,length(noiseVec));
    cubicDep = zeros(2,length(noiseVec));
    sinusoidalDep = zeros(2,length(noiseVec));
    hiFreqSinDep = zeros(2,length(noiseVec));
    fourthRootDep = zeros(2,length(noiseVec));
    circleDep = zeros(2,length(noiseVec));
    stepDep = zeros(2,length(noiseVec));
    
    noiseVecIdx = 1:length(noiseVec);
    
    % collect all the data for each of the scanincrs we tested and store in
    % a matrix, for each dependency type, we plot the difference between
    % the minimum and the maximum
    for ii=1:length(scanincrsToTest)
        scanincrVal = scanincrsToTest(ii);
        % TODO: INDEX PROPERLY -- currently assumes contiguous noiseVec
        linearData     = squeeze(powerCurve(ii,1,noiseVecIdx));
        quadraticData  = squeeze(powerCurve(ii,2,noiseVecIdx));
        cubicData      = squeeze(powerCurve(ii,3,noiseVecIdx));
        sinusoidalData = squeeze(powerCurve(ii,4,noiseVecIdx));
        hiFreqSinData  = squeeze(powerCurve(ii,5,noiseVecIdx));
        fourthRootData = squeeze(powerCurve(ii,6,noiseVecIdx));
        circleData     = squeeze(powerCurve(ii,7,noiseVecIdx));
        stepData       = squeeze(powerCurve(ii,8,noiseVecIdx));
        
        % seed the data
        if(ii==1)
            for jj=1:2
                linearDep(jj,:) = linearData';
                quadraticDep(jj,:) = quadraticData';
                cubicDep(jj,:) = cubicData';
                fourthRootDep(jj,:) = fourthRootData';
                circleDep(jj,:) = circleData';
                stepDep(jj,:) = stepData';
            end
        end
        if(scanincrVal<0.1138)
            for jj=1:2
                sinusoidalDep(jj,:) = sinusoidalData;
            end
        end
        if(scanincrVal<0.029)
            for jj=1:2
                hiFreqSinDep(jj,:) = hiFreqSinData;
            end
        end

        for jj=1:length(noiseVec)
            % process linear
            if(linearData(jj)<linearDep(1,jj))
                linearDep(1,jj) = linearData(jj);
            end
            if(linearData(jj)>linearDep(2,jj))
                linearDep(2,jj) = linearData(jj);
            end

            % process quadratic
            if(quadraticData(jj)<quadraticDep(1,jj))
                quadraticDep(1,jj) = quadraticData(jj);
            end
            if(quadraticData(jj)>quadraticDep(2,jj))
                quadraticDep(2,jj) = quadraticData(jj);
            end

            % process cubic
            if(cubicData(jj)<cubicDep(1,jj))
                cubicDep(1,jj) = cubicData(jj);
            end
            if(cubicData(jj)>cubicDep(2,jj))
                cubicDep(2,jj) = cubicData(jj);
            end

            % process sinusoidal
            if(scanincrVal<0.1138)
                if(sinusoidalData(jj)<sinusoidalDep(1,jj))
                    sinusoidalDep(1,jj) = sinusoidalData(jj);
                end
                if(sinusoidalData(jj)>sinusoidalDep(2,jj))
                    sinusoidalDep(2,jj) = sinusoidalData(jj);
                end
            end

            % process hi-freq sine
            if(scanincrVal<0.029)
                if(hiFreqSinData(jj)<hiFreqSinDep(1,jj))
                    hiFreqSinDep(1,jj) = hiFreqSinData(jj);
                end
                if(hiFreqSinData(jj)>hiFreqSinDep(2,jj))
                    hiFreqSinDep(2,jj) = hiFreqSinData(jj);
                end
            end

            % process fourth-root
            if(fourthRootData(jj)<fourthRootDep(1,jj))
                fourthRootDep(1,jj) = fourthRootData(jj);
            end
            if(fourthRootData(jj)>fourthRootDep(2,jj))
                fourthRootDep(2,jj) = fourthRootData(jj);
            end

            % process circular
            if(circleData(jj)<circleDep(1,jj))
                circleDep(1,jj) = circleData(jj);
            end
            if(circleData(jj)>circleDep(2,jj))
                circleDep(2,jj) = circleData(jj);
            end

            % process step
            if(stepData(jj)<stepDep(1,jj))
                stepDep(1,jj) = stepData(jj);
            end
            if(stepData(jj)>stepDep(2,jj))
                stepDep(2,jj) = stepData(jj);
            end

        end
        numLoopsRun = numLoopsRun + 1;
    end
    
    linearDepCell{MIdx}     = linearDep;
    quadraticDepCell{MIdx}  = quadraticDep;
    cubicDepCell{MIdx}      = cubicDep;
    sinusoidalDepCell{MIdx} = sinusoidalDep;
    hiFreqSinDepCell{MIdx}  = hiFreqSinDep;
    fourthRootDepCell{MIdx} = fourthRootDep;
    circleDepCell{MIdx}     = circleDep;
    stepDepCell{MIdx}       = stepDep;
end

linearDepToPlot = zeros(2,length(noiseVec));
quadraticDepToPlot = zeros(2,length(noiseVec));
cubicDepToPlot = zeros(2,length(noiseVec));
sinusoidalDepToPlot = zeros(2,length(noiseVec));
hiFreqSinDepToPlot = zeros(2,length(noiseVec));
fourthRootDepToPlot = zeros(2,length(noiseVec));
circleDepToPlot = zeros(2,length(noiseVec));
stepDepToPlot = zeros(2,length(noiseVec));

for ii=1:length(noiseVec)
    linearDepToPlot(1,ii) = linearDepCell{1}(1,ii);
    linearDepToPlot(2,ii) = linearDepCell{1}(2,ii);
    quadraticDepToPlot(1,ii) = quadraticDepCell{1}(1,ii);
    quadraticDepToPlot(2,ii) = quadraticDepCell{1}(2,ii);
    cubicDepToPlot(1,ii) = cubicDepCell{1}(1,ii);
    cubicDepToPlot(2,ii) = cubicDepCell{1}(2,ii);
    sinusoidalDepToPlot(1,ii) = sinusoidalDepCell{1}(1,ii);
    sinusoidalDepToPlot(2,ii) = sinusoidalDepCell{1}(2,ii);
    hiFreqSinDepToPlot(1,ii) = hiFreqSinDepCell{1}(1,ii);
    hiFreqSinDepToPlot(2,ii) = hiFreqSinDepCell{1}(2,ii);
    fourthRootDepToPlot(1,ii) = fourthRootDepCell{1}(1,ii);
    fourthRootDepToPlot(2,ii) = fourthRootDepCell{1}(2,ii);
    circleDepToPlot(1,ii) = circleDepCell{1}(1,ii);
    circleDepToPlot(2,ii) = circleDepCell{1}(2,ii);
    stepDepToPlot(1,ii) = stepDepCell{1}(1,ii);
    stepDepToPlot(2,ii) = stepDepCell{1}(2,ii);
    
    for jj=2:length(MVecToPlot)
        rangeCur  = linearDepToPlot(2,ii)-linearDepToPlot(1,ii);
        rangePrev = linearDepCell{jj-1}(2,ii)-linearDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            linearDepToPlot(1,ii) = linearDepCell{jj}(1,ii);
            linearDepToPlot(2,ii) = linearDepCell{jj}(2,ii);
        end
        
        rangeCur  = quadraticDepToPlot(2,ii)-quadraticDepToPlot(1,ii);
        rangePrev = quadraticDepCell{jj-1}(2,ii)-quadraticDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            quadraticDepToPlot(1,ii) = quadraticDepCell{jj}(1,ii);
            quadraticDepToPlot(2,ii) = quadraticDepCell{jj}(2,ii);
        end
        
        rangeCur  = cubicDepToPlot(2,ii)-cubicDepToPlot(1,ii);
        rangePrev = cubicDepCell{jj-1}(2,ii)-cubicDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            cubicDepToPlot(1,ii) = cubicDepCell{jj}(1,ii);
            cubicDepToPlot(2,ii) = cubicDepCell{jj}(2,ii);
        end
        
        rangeCur  = sinusoidalDepToPlot(2,ii)-sinusoidalDepToPlot(1,ii);
        rangePrev = sinusoidalDepCell{jj-1}(2,ii)-sinusoidalDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            sinusoidalDepToPlot(1,ii) = sinusoidalDepCell{jj}(1,ii);
            sinusoidalDepToPlot(2,ii) = sinusoidalDepCell{jj}(2,ii);
        end
        
        rangeCur  = hiFreqSinDepToPlot(2,ii)-hiFreqSinDepToPlot(1,ii);
        rangePrev = hiFreqSinDepCell{jj-1}(2,ii)-hiFreqSinDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            hiFreqSinDepToPlot(1,ii) = hiFreqSinDepCell{jj}(1,ii);
            hiFreqSinDepToPlot(2,ii) = hiFreqSinDepCell{jj}(2,ii);
        end
        
        rangeCur  = fourthRootDepToPlot(2,ii)-fourthRootDepToPlot(1,ii);
        rangePrev = fourthRootDepCell{jj-1}(2,ii)-fourthRootDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            fourthRootDepToPlot(1,ii) = fourthRootDepCell{jj}(1,ii);
            fourthRootDepToPlot(2,ii) = fourthRootDepCell{jj}(2,ii);
        end
        
        rangeCur  = circleDepToPlot(2,ii)-circleDepToPlot(1,ii);
        rangePrev = circleDepCell{jj-1}(2,ii)-circleDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            circleDepToPlot(1,ii) = circleDepCell{jj}(1,ii);
            circleDepToPlot(2,ii) = circleDepCell{jj}(2,ii);
        end
        
        rangeCur  = stepDepToPlot(2,ii)-stepDepToPlot(1,ii);
        rangePrev = stepDepCell{jj-1}(2,ii)-stepDepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            stepDepToPlot(1,ii) = stepDepCell{jj}(1,ii);
            stepDepToPlot(2,ii) = stepDepCell{jj}(2,ii);
        end
    end
end

noiseVecToPlot = noiseVec/10;

figure;
h = subplot(2,4,1);
plot(noiseVecToPlot,linearDepToPlot(2,:)-linearDepToPlot(1,:));
grid on; 
title('Linear','FontSize',20);
h.FontSize = 20;

h = subplot(2,4,2);
plot(noiseVecToPlot,quadraticDepToPlot(2,:)-quadraticDepToPlot(1,:));
grid on; 
title('Quadratic','FontSize',20);
h.FontSize = 20;

h = subplot(2,4,3);
plot(noiseVecToPlot,cubicDepToPlot(2,:)-cubicDepToPlot(1,:));
grid on; 
title('Cubic','FontSize',20);
h.FontSize = 20;

h = subplot(2,4,4);
plot(noiseVecToPlot,sinusoidalDepToPlot(2,:)-sinusoidalDepToPlot(1,:));
grid on; 
title('Sinusoidal','FontSize',20);
h.FontSize = 20;

h = subplot(2,4,5);
plot(noiseVecToPlot,hiFreqSinDepToPlot(2,:)-hiFreqSinDepToPlot(1,:));
grid on; 
title('Hi-Freq Sin','FontSize',20); 
h.FontSize = 20;

h = subplot(2,4,6);
plot(noiseVecToPlot,fourthRootDepToPlot(2,:)-fourthRootDepToPlot(1,:));
grid on; 
title('Fourth-Root','FontSize',20);
h.FontSize = 20;

h = subplot(2,4,7);
plot(noiseVecToPlot,circleDepToPlot(2,:)-circleDepToPlot(1,:));
grid on; 
title('Circular','FontSize',20);
h.FontSize = 20;

h = subplot(2,4,8);
plot(noiseVecToPlot,stepDepToPlot(2,:)-stepDepToPlot(1,:));
grid on; 
title('Step','FontSize',20);
h.FontSize = 20;

[~,h] = suplabel('Noise','x');
set(h,'FontSize',20);
[~,h] = suplabel('$$max[\Delta \Gamma(\widehat{CIM})]$$','y');
h.FontSize = 20;
h.Interpreter = 'Latex';

end