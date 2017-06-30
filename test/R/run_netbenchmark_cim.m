%% Script which processes the netbenchmark data
clear;
clc;

% setup and start the parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 6;
saveProfile(myCluster);
p = gcp

dataRepo = '/home/kiran/ownCloud/PhD/sim_results/netbenchmark/data/';
dataSourcesToProcess = {'syntren300','rogers1000','syntren1000','gnw1565','gnw2000'};
numDatasets = 5;

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
        
        fOut = fullfile(dataRepo,sprintf('%s_%d_output.mat',dataSource,ii));
        save(fOut,'R');
    end
end