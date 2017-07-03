%% Script which processes the netbenchmark data
clear;
clc;

% setup and start the parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 6;
saveProfile(myCluster);
p = gcp

dataRepo = '/home/kiran/data/netbenchmark/inputs/';
dataSourcesToProcess = {'syntren300','rogers1000','syntren1000','gnw1565','gnw2000'};
dataOutputFolder = '/home/kiran/data/netbenchmark/matlab_outputs';
numDatasets = 150;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for dataSourceIdx=1:length(dataSourcesToProcess)
    dataSource = dataSourcesToProcess{dataSourceIdx};
    for ii=1:numDatasets
        dispstat(sprintf('Computing %s[%d/%d]',dataSource, ii, numDatasets),'keepthis', 'timestamp');
        
        fIn = fullfile(dataRepo,sprintf('%s_%d.mat',dataSource,ii));
        load(fIn);
        
        X = data;
        R = paircim_v4_cc_mexoffload( X );
        
        fOut = fullfile(dataOutputFolder,sprintf('%s_%d_output.mat',dataSource,ii));
        save(fOut,'R');
    end
end

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