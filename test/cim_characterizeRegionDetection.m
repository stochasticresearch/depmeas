%% Characterize the region detection performance of CIM

clear;
clc;
close all;
dbstop if error;

rng(1234);

numMCSims = 500;
MVec = [50 100 200:100:1000];
noiseVec = 0:30;
% simulate different change rates
cornerPts = 0.1:0.1:0.9;
minScanIncr = 0.015625;

rdData_up = zeros(length(MVec),length(noiseVec),length(cornerPts),numMCSims);
rdData_down = zeros(length(MVec),length(noiseVec),length(cornerPts),numMCSims);

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for mIdx=1:length(MVec)
    M = MVec(mIdx);
    for noiseIdx=1:length(noiseVec)
        l = noiseVec(noiseIdx);
        dispstat(sprintf('Computing for M=%d noise=%d', M, l),'keepthis', 'timestamp');
        for cornerPtIdx=1:length(cornerPts)
            cornerPt = cornerPts(cornerPtIdx);
            parfor ii=1:numMCSims
                x = rand(M,1);
                y1 = 4*(x-cornerPt).^2 + noise*(l/num_noise)*randn(M,1);
                y2 = -4*(x-cornerPt).^2 + noise*(l/num_noise)*randn(M,1);
                
                [~,rect1] = cim(x,y1,minScanIncr);
                [~,rect2] = cim(x,y1,minScanIncr);
                
                % parse where the region detection found the division
                % marker for the down parabola
                
                if(rect1(4,1)==1)
                    % assume we are u/v orientation
                    r1 = rect1(2,1);
                else
                    % assume we are v/u orientation
                    r1 = rect1(4,1);
                end
                
                % parse where the region detection found the division
                % marker for the up parabola
                if(rect2(4,1)==1)
                    % assume we are u/v orientation
                    r2 = rect2(2,1);
                else
                    % assume we are v/u orientation
                    r2 = rect2(4,1);
                end
                fprintf('r1=%0.02f r2=%0.02f\n', r1, r2);
                rdData_up(mIdx,noiseIdx,cornerPtIdx,ii) = r1;
                rdData_down(mIdx,noiseIdx,cornerPtIdx,ii) = r2;
            end
        end
    end
    % store the results
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_regionDetection_%d.mat',M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_regionDetection_%d.mat',M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_regionDetection_%d.mat',M));
    end
end

%% Plot the results
clear;
clc;

MVecOfInterest = [100 200 300 400 500];
noiseOfInterest = [0 10 20];
cornerPtsOfInterest = [0.2 0.5 0.7];
fontSizeVal = 20;

M_complete = 700;  % change this tot he latest run that is complete ...
                   % we unfortuantely did this badly ... the way we saved
                   % the data
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_regionDetection_%d.mat',M_complete));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_regionDetection_%d.mat',M_complete));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_regionDetection_%d.mat',M_complete));
end

    
p = numSubplots(length(noiseOfInterest)*length(cornerPtsOfInterest));
plotIdx = 1;
for cornerPtIdx=1:length(cornerPtsOfInterest)
    % find the correct index into the rdData matrix
    cornerPt = cornerPtsOfInterest(cornerPtIdx);
    cornerPt_ii = find(cornerPts==cornerPt);
    for noiseIdx=1:length(noiseOfInterest)
        noise = noiseOfInterest(noiseIdx);
        noise_ii = find(noiseVec==noise);
        
        boxData = zeros(numMCSims,length(MVecOfInterest));
        labelCell = cell(1,length(MVecOfInterest));
        ii = 1;
        for M=MVecOfInterest
            MIdx = find(M==MVec);
            upData = squeeze(rdData_up(MIdx,noise_ii,cornerPt_ii,:));
            boxData(:,ii) = upData;
            labelCell{ii} = sprintf('M=%d',M);
            ii = ii + 1;
        end
        subplot(p(1),p(2),plotIdx);
        if(plotIdx<=3)
            bh=boxplot(boxData,'Labels',labelCell);
            if(plotIdx==2)
                title(sprintf('Noise=%d',noise/10),'FontSize',fontSizeVal);
            else
                title(sprintf('Noise=%d',noise/10),'FontSize',fontSizeVal);
            end
        else
            bh=boxplot(boxData,'Labels',{'','','','',''});
        end

        ylim([0 1])
        set(bh,'LineWidth',2.5);
        plotIdx = plotIdx + 1;
        hh = gca;
        hh.FontSize = fontSizeVal;
        hold on;
        plot(xlim,[cornerPt cornerPt],'g:','LineWidth',2);
    end
end
    
