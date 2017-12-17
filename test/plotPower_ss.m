function [] = plotPower_ss(sampleSizeAnalysisVec, labels, noiseVec, style)

numDepTypes = 8;
M_inlet = 200;
% generate the inlet data
inletX = linspace(0,1,M_inlet);
inletT = linspace(0,2*pi,M_inlet);
inletData = zeros(numDepTypes,M_inlet);
inletData(1,:) = inletX;
inletData(2,:) = 4*(inletX-.5).^2;
inletData(3,:) = 128*(inletX-1/3).^3-48*(inletX-1/3).^3-12*(inletX-1/3);
inletData(4,:) = sin(4*pi*inletX);
inletData(5,:) = sin(16*pi*inletX);
inletData(6,:) = inletX.^(1/4);
inletData(8,:) = (inletX > 0.5);
        
if(style==1)
    % the original Simon/Tibshirani style of power-plots
    plotStyle = {'o-.', '+-.', 'd-.', 'v-.', 's-.', 'p-.'};
    captionText = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'};
    % inlet plot configuration
    
    inset_bufX = 0.0; inset_bufY = 0.26;
    inset_width = 0.1; inset_height = 0.08;
    
    sampleSizeAnalysisVec(sampleSizeAnalysisVec==0)=NaN;

    % we break it up into 2 figures, 2x2 subplots on each figure for
    % readability
    captionTextIdx = 1;
    for figNum=[1,2]
        figure(figNum);
        for subplotNum=1:4
            depTypeIdx = (figNum-1)*4+subplotNum;
            hhCell = cell(1,length(labels));
            for ii=1:length(labels)
                h = subplot(2,2,subplotNum);
                hh = plot(noiseVec, squeeze(sampleSizeAnalysisVec(ii,depTypeIdx,1:length(noiseVec))), plotStyle{ii});
                hhCell{ii} = hh;
                hold on;
            end
            xlabel(captionText{captionTextIdx}, 'FontSize', 20); 
            captionTextIdx = captionTextIdx + 1;
            grid on;
            h.FontSize = 20; 
            loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
            ax = axes('Position',loc_inset);
            if(depTypeIdx~=7)
                plot(inletX,inletData(depTypeIdx,:), 'k', 'LineWidth', 2);
                ax.XLim = [min(inletX) max(inletX)];
                ax.YLim = [min(inletData(depTypeIdx,:)) max(inletData(depTypeIdx,:))];
            else
                plot(cos(inletT),sin(inletT), 'k', 'LineWidth', 2);
                ax.XLim = [min(cos(inletT)) max(cos(inletT))];
                ax.YLim = [min(sin(inletT)) max(sin(inletT))];
            end
            ax.Box = 'on'; ax.XTick = []; ax.YTick = [];
            hhCell{1}.LineWidth = 5; 
            for ii=2:length(labels)
                hhCell{ii}.LineWidth = 1.5; 
            end
        end
        if(figNum==2)
            legend(h,labels);
        end
        [~,h] = suplabel('Noise','x');
        set(h,'FontSize',20);
        [~,h] = suplabel('$$M_{min}$$','y');
        h.FontSize = 20;
        h.Interpreter = 'Latex';
    end
elseif(style==2)
    numRelationshipsPerFig = 4;
    numVertSubplot = 4;
    cmap = helper_power_colormap();
    
    width = 50; height = width/5;
    font_size = 15;
    
    for figNum=[1,2]
        figure('paperpositionmode', 'auto', 'units', 'centimeters', 'position', [0 0 width height])
        
        for subplotNum=1:numRelationshipsPerFig
            depTypeIdx = (figNum-1)*4+subplotNum;
            % plot the relationship type
            subplot(numVertSubplot,numRelationshipsPerFig,subplotNum);
            if(depTypeIdx~=7)
                plot(inletX,inletData(depTypeIdx,:), '.k', 'markerfacecolor', 'k', 'markersize', 4);
            else
                plot(cos(inletT),sin(inletT), '.k', 'markerfacecolor', 'k', 'markersize', 4);
            end
            set(gca, 'box', 'on', 'xtick', [], 'ytick', [])
            axis image
            axis square
            axis off
            
            % plot the power
            % create teh subplotVec
            spVec = [];
            cVal = subplotNum;
            for zz=2:numVertSubplot
                cVal = cVal+numRelationshipsPerFig;
                spVec = [spVec cVal];
            end
            subplot(numVertSubplot,4,spVec)
            pm = squeeze(sampleSizeAnalysisVec(:,depTypeIdx,1:length(noiseVec)));
            pm(pm==0)=inf;
            imagesc(pm,'AlphaData',pm);
            
            if(subplotNum==1)
                set(gca, 'ytick', 1:length(labels), 'yticklabel', labels, 'FontSize', font_size)
            else
                set(gca, 'ytick', 1:length(labels), 'yticklabel', [], 'FontSize', font_size)
            end
            
            if(subplotNum==numRelationshipsPerFig && figNum==2)
                h = colorbar();
                set(h,'fontsize',font_size, 'position', [0.9256    0.2303    0.0125    0.2602], 'linewidth', 0.5);
                ylabel(h,'# Samples')
            end
            
        end
        if(figNum==2)
            [~,h3]=suplabel('Noise','x',[.0975 .12 .84 .84]);
            set(h3,'FontSize',font_size)
        end
        colormap(cmap)
    end
    
else
    error('Current style not supported!!');
end

end
