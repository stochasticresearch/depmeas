%% Plot the netbenchmark combined data
clear;
clc;
close all;

datasetsToPlot = {'syntren300'};
if(ispc)
    combinedDataRepo = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\netbenchmark';
elseif(ismac)
    combinedDataRepo = '/Users/Kiran/ownCloud/PhD/sim_results/netbenchmark';
else
    combinedDataRepo = '/home/kiran/ownCloud/PhD/sim_results/netbenchmark';
end

estimatorsToPlot = {'"CIM"', '"Kendall"', '"KNN1"', '"KNN6"', '"KNN20"', '"vME"', '"AP"', '"rand"'};
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
fontSizeVal = 20;
for datasetIdx=1:length(datasetsToPlot)
    figure;

    ds = datasetsToPlot{datasetIdx};
    combinedDataFile = fullfile(combinedDataRepo,sprintf('%s_combined.csv',ds));
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
    
    % box plot the results for each column
    h = boxplot(x,'Labels',y,'Notch','on');
    title(sprintf('%s',ds), 'FontSize', fontSizeVal);
    grid on;
    ylabel('max[AUPR_{20}]', 'FontSize', fontSizeVal);
    hh = gca;
    hh.FontSize = fontSizeVal;
    set(h,'LineWidth',3);
end