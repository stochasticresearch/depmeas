%% Plot the netbenchmark combined data
clear;
clc;
close all;

datasetsToPlot = {'syntren300'};
totalEdgesSyntren300=468;
globalNoiseVec = [0,10,20,30,40,50];
localNoiseVec = [5,10,20,30,40,50];

% globalNoiseVec = [10];
% localNoiseVec = [5];

if(ispc)
    combinedDataRepo = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\netbenchmark';
elseif(ismac)
    combinedDataRepo = '/Users/kiran/ownCloud/PhD/sim_results/netbenchmark';
else
    combinedDataRepo = '/home/kiran/ownCloud/PhD/sim_results/netbenchmark';
end

%estimatorsToPlot = {'"CIM"', '"Kendall"', '"KNN1"', '"KNN6"', '"KNN20"', '"vME"', '"AP"', '"rand"'};
estimatorsToPlot = {'"CIM"', '"KNN1"', '"KNN6"', '"KNN20"', '"vME"', '"AP"'};
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
fontSizeVal = 22;
piVec = [];
for datasetIdx=1:length(datasetsToPlot)
    for gn=globalNoiseVec
        for ln=localNoiseVec
            % if the file exists -- process it
            ds = datasetsToPlot{datasetIdx};
            combinedDataFile = fullfile(combinedDataRepo,sprintf('%s_combined_gn_%d_ln_%d.csv',ds,gn,ln));
            if(exist(combinedDataFile,'file'))
                A = importdata(combinedDataFile);

                % select the correct columns from the data, based on which estimators
                % are desired to be plotted
                idxs = zeros(1,length(estimatorsToPlot));
                for ii=1:length(estimatorsToPlot)
                    fIdx = find(cellfun(cellfind(estimatorsToPlot{ii}),A.colheaders));
                    idxs(ii) = fIdx;
                end
                x = A.data(:,idxs);
                y = strrep(A.colheaders(idxs),'"','');
                
                if(gn==10 && ln==5)  % plot best case scenario :0
                    figure;
                    % box plot the results for each column
                    h = boxplot(x,'Labels',y,'Notch','on');
%                     title(sprintf('%s - [%d/%d]',ds,gn,ln), 'FontSize', fontSizeVal);
                    title(sprintf('%s',ds), 'FontSize', fontSizeVal);
                    grid on;
                    ylabel('max[AUPR_{20}]', 'FontSize', fontSizeVal);
                    hh = gca;
                    hh.FontSize = fontSizeVal;
                    set(h,'LineWidth',3);
                    set(gca, 'YScale', 'log');
                end
                
                medX = median(x);
                cimMed = medX(1); nextBest = max(medX(2:end));
                percentageImprovement = (cimMed-nextBest)/nextBest;
                piVec = [piVec percentageImprovement];
                
                fprintf('**********\n');
                g=sprintf('%0.04f ', medX);
                fprintf('gn=%d ln=%d >> %s << %0.04f\n', gn,ln,g,percentageImprovement);
                fprintf('**********\n');
            end
        end
    end
end

% print out some statistics of the piVec
fprintf('----------------------------------------------------------\n');
fprintf('minImprovement=%0.02f maxImprovement=%.02f meanImprovement=%0.02f\n',...
        min(piVec)*100,max(piVec)*100,mean(piVec)*100);
fprintf('minEdges=%d maxEdges=%d meanEdges=%d\n',...
        floor(min(piVec)*totalEdgesSyntren300),...
        floor(max(piVec)*totalEdgesSyntren300),...
        floor(mean(piVec)*totalEdgesSyntren300));
fprintf('----------------------------------------------------------\n');