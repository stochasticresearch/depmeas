function [] = plotAlgoAlphaSensitivity_withinM(MVecToPlot,num_noise_test_min,num_noise_test_max)

fnameStr = 'cim_v2_cc';
CIMVECIDX = 3;
noiseVec = num_noise_test_min:num_noise_test_max;

linearDepCell = cell(1,length(MVecToPlot));
quadraticDepCell = cell(1,length(MVecToPlot));
cubicDepCell = cell(1,length(MVecToPlot));
sinusoidalDepCell = cell(1,length(MVecToPlot));
hiFreqSinDepCell = cell(1,length(MVecToPlot));
fourthRootDepCell = cell(1,length(MVecToPlot));
circleDepCell = cell(1,length(MVecToPlot));
stepDepCell = cell(1,length(MVecToPlot));
indepCell = cell(1,length(MVecToPlot));

numLoopsRun = 1;
for MIdx=1:length(MVecToPlot)
    M = MVecToPlot(MIdx);
    % load the data
    if(ispc)
        load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_algoSensitivityAlpha_M_%d.mat', fnameStr, M));
    elseif(ismac)
        load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivityAlpha_M_%d.mat', fnameStr, M));
    else
        load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivityAlpha_M_%d.mat', fnameStr, M));
    end
    
    linearDep = zeros(2,length(noiseVec));      % 1 - min
                                                % 2 - max
    quadraticDep = zeros(2,length(noiseVec));
    cubicDep = zeros(2,length(noiseVec));
    sinusoidalDep = zeros(2,length(noiseVec));
    hiFreqSinDep = zeros(2,length(noiseVec));
    fourthRootDep = zeros(2,length(noiseVec));
    circleDep = zeros(2,length(noiseVec));
    stepDep = zeros(2,length(noiseVec));
    indep = zeros(2,length(noiseVec));
    
    noiseVecIdx = 1:length(noiseVec);
   
%     alphasToTest = [0.2];
    
    % collect all the data for each of the scanincrs we tested and store in
    % a matrix, for each dependency type, we plot the difference between
    % the minimum and the maximum
    for ii=1:length(alphasToTest)
        rawData = algoSensitivityData{ii};  % algoSensitivityData is loaded from the file above
        linearData = rawData.linearDep(CIMVECIDX,noiseVecIdx);
        quadraticData = rawData.quadraticDep(CIMVECIDX,noiseVecIdx);
        cubicData = rawData.cubicDep(CIMVECIDX,noiseVecIdx);
        sinusoidalData = rawData.sinusoidalDep(CIMVECIDX,noiseVecIdx);
        hiFreqSinData = rawData.hiFreqSinDep(CIMVECIDX,noiseVecIdx);
        fourthRootData = rawData.fourthRootDep(CIMVECIDX,noiseVecIdx);
        circleData = rawData.circleDep(CIMVECIDX,noiseVecIdx);
        stepData = rawData.stepDep(CIMVECIDX,noiseVecIdx);
        indepData = rawData.indep(CIMVECIDX,noiseVecIdx);

        % seed the data
        if(ii==1)
            for jj=1:2
                linearDep(jj,noiseVecIdx) = linearData;
                quadraticDep(jj,noiseVecIdx) = quadraticData;
                cubicDep(jj,noiseVecIdx) = cubicData;
                sinusoidalDep(jj,noiseVecIdx) = sinusoidalData;
                hiFreqSinDep(jj,noiseVecIdx) = hiFreqSinData;
                fourthRootDep(jj,noiseVecIdx) = fourthRootData;
                circleDep(jj,noiseVecIdx) = circleData;
                cubicDep(jj,noiseVecIdx) = cubicData;
                stepDep(jj,noiseVecIdx) = stepData;
                indep(jj,noiseVecIdx) = indepData;                
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
            if(sinusoidalData(jj)<sinusoidalDep(1,jj))
                sinusoidalDep(1,jj) = sinusoidalData(jj);
            end
            if(sinusoidalData(jj)>sinusoidalDep(2,jj))
                sinusoidalDep(2,jj) = sinusoidalData(jj);
            end

            % process hi-freq sine
            if(hiFreqSinData(jj)<hiFreqSinDep(1,jj))
                hiFreqSinDep(1,jj) = hiFreqSinData(jj);
            end
            if(hiFreqSinData(jj)>hiFreqSinDep(2,jj))
                hiFreqSinDep(2,jj) = hiFreqSinData(jj);
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

            % process indep
            if(indepData(jj)<indep(1,jj))
                indep(1,jj) = indepData(jj);
            end
            if(indepData(jj)>indep(2,jj))
                indep(2,jj) = indepData(jj);
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
    indepCell{MIdx}         = indep;    
end

linearDepToPlot = zeros(2,length(noiseVec));
quadraticDepToPlot = zeros(2,length(noiseVec));
cubicDepToPlot = zeros(2,length(noiseVec));
sinusoidalDepToPlot = zeros(2,length(noiseVec));
hiFreqSinDepToPlot = zeros(2,length(noiseVec));
fourthRootDepToPlot = zeros(2,length(noiseVec));
circleDepToPlot = zeros(2,length(noiseVec));
stepDepToPlot = zeros(2,length(noiseVec));
indepToPlot = zeros(2,length(noiseVec));

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
    indepToPlot(1,ii) = indepCell{1}(1,ii);
    indepToPlot(2,ii) = indepCell{1}(2,ii);
    
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
        
        rangeCur  = indepToPlot(2,ii)-indepToPlot(1,ii);
        rangePrev = indepCell{jj-1}(2,ii)-indepCell{jj-1}(1,ii);
        if(rangeCur>rangePrev)
            indepToPlot(1,ii) = indepCell{jj}(1,ii);
            indepToPlot(2,ii) = indepCell{jj}(2,ii);
        end
    end
end

% generate the inlet data
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

noiseVecToPlot = noiseVec/10;
lineWidthVal = 2.5;

y1 = linearDepToPlot(2,:)-linearDepToPlot(1,:);
y2 = quadraticDepToPlot(2,:)-quadraticDepToPlot(1,:);
y3 = cubicDepToPlot(2,:)-cubicDepToPlot(1,:);
y4 = sinusoidalDepToPlot(2,:)-sinusoidalDepToPlot(1,:);
y5 = hiFreqSinDepToPlot(2,:)-hiFreqSinDepToPlot(1,:);
y6 = fourthRootDepToPlot(2,:)-fourthRootDepToPlot(1,:);
y7 = circleDepToPlot(2,:)-circleDepToPlot(1,:);
y8 = stepDepToPlot(2,:)-stepDepToPlot(1,:);

minX = min(noiseVecToPlot);
maxX = max(noiseVecToPlot);
minY = min([y1 y2 y3 y4 y5 y6 y7 y8]);
maxY = max([y1 y2 y3 y4 y5 y6 y7 y8]);

if(minY==maxY)
    maxY = minY + eps;
end

inset_bufX = 0; inset_bufY = 0.28;
inset_width = 0.06; inset_height = 0.06;

% figure;
width=30;
height=width/3.2;
figure('paperpositionmode', 'auto', 'units', 'centimeters', 'position', [0 0 width height])
h = subplot(2,4,1);
plot(noiseVecToPlot,y1,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Linear','FontSize',20);
h.FontSize = 20;
h.YTick = [0, 0.02, 0.04, 0.06];
inletIdx = 1;
loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
ax.XLim = [min(inletX) max(inletX)];
ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

h = subplot(2,4,2);
plot(noiseVecToPlot,y2,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Quadratic','FontSize',20); 
h.YTick = [0, 0.02, 0.04, 0.06];
set(h,'yticklabel',[])
h.FontSize = 20;
inletIdx = 2;
loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
ax.XLim = [min(inletX) max(inletX)];
ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

h = subplot(2,4,3);
plot(noiseVecToPlot,y3,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Cubic','FontSize',20);
% legend({'Correction','No Correction'});
h.YTick = [0, 0.02, 0.04, 0.06];
set(h,'yticklabel',[])
h.FontSize = 20;
inletIdx = 3;
loc_inset = [h.Position(1)+inset_bufX+0.095 h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
ax.XLim = [min(inletX) max(inletX)];
ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

h = subplot(2,4,4);
plot(noiseVecToPlot,y4,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Sinusoidal','FontSize',20); 
% legend({'Correction','No Correction'});
h.YTick = [0, 0.02, 0.04, 0.06];
set(h,'yticklabel',[])
h.FontSize = 20;
inletIdx = 4;
loc_inset = [h.Position(1)+inset_bufX+0.095 h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
ax.XLim = [min(inletX) max(inletX)];
ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

h = subplot(2,4,5);
plot(noiseVecToPlot,y5,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Hi-Freq Sin','FontSize',20); 
% legend({'Correction','No Correction'});
set(h,'xticklabel',[])
h.YTick = [0, 0.02, 0.04, 0.06];
h.FontSize = 20;
inletIdx = 5;
loc_inset = [h.Position(1)+inset_bufX+0.095 h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
ax.XLim = [min(inletX) max(inletX)];
ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

h = subplot(2,4,6);
plot(noiseVecToPlot,y6,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Fourth-Root','FontSize',20);
h.YTick = [0, 0.02, 0.04, 0.06];
set(h,'xticklabel',[],'yticklabel',[])
h.FontSize = 20;
inletIdx = 6;
loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
ax.XLim = [min(inletX) max(inletX)];
ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

h = subplot(2,4,7);
plot(noiseVecToPlot,y7,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Circular','FontSize',20);
h.YTick = [0, 0.02, 0.04, 0.06];
set(h,'xticklabel',[],'yticklabel',[])
h.FontSize = 20;
inletIdx = 7;
loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(cos(inletT),sin(inletT), 'k', 'LineWidth', 2);
ax.XLim = [min(cos(inletT)) max(cos(inletT))];
ax.YLim = [min(sin(inletT)) max(sin(inletT))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

h = subplot(2,4,8);
plot(noiseVecToPlot,y8,'LineWidth',lineWidthVal); axis([minX maxX minY maxY]);
grid on; 
% title('Step','FontSize',20);
h.YTick = [0, 0.02, 0.04, 0.06];
set(h,'xticklabel',[],'yticklabel',[])
h.FontSize = 20;
inletIdx = 8;
loc_inset = [h.Position(1)+inset_bufX+0.095 h.Position(2)+inset_bufY inset_width inset_height];
ax = axes('Position',loc_inset);
plot(inletX,inletData(inletIdx,:), 'k', 'LineWidth', 2);
ax.XLim = [min(inletX) max(inletX)];
ax.YLim = [min(inletData(inletIdx,:)) max(inletData(inletIdx,:))];
ax.Box = 'on'; ax.XTick = []; ax.YTick = [];

[~,h] = suplabel('Noise','x',[.08 .15 .84 .84]);
set(h,'FontSize',20);
[~,h] = suplabel('$$max[\Delta \widehat{CIM}]$$','y');
h.FontSize = 20;
h.Interpreter = 'Latex';

% subplot(3,3,9);
% plot(noiseVec,indepToPlot(2,:)-indepToPlot(1,:));
% grid on; xlabel('Noise'); title('Independence'); ylabel('max[$$\Delta \widehat{CIM}]$$','interpreter','Latex');

end