%% Plot the netbenchmark combined data
clear;
clc;
close all;

datasetsToPlot = {'syntren300'};
if(ispc)
    error('fill in!')
elseif(ismac)
    combinedDataRepo = '/Users/Kiran/data/netbenchmark/combined_outputs';
else
    combinedDataRepo = '/home/kiran/data/netbenchmark/combined_outputs';
end

for datasetIdx=1:length(datasetsToPlot)
    figure;

    ds = datasetsToPlot{datasetIdx};
    combinedDataFile = fullfile(combinedDataRepo,sprintf('%s_combined.csv',ds));
    A = importdata(combinedDataFile);
    
    % box plot the results for each column
    boxplot(A.data,'Labels',A.colheaders);
    title(sprintf('%s',ds));
    grid on;
    ylabel('max(AUPR)')

end