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
        [R_cim,R_knn1,R_knn6,R_knn20,R_vmeMI,R_apMI] = pair_all_cc_mex( X );
        
        fOut = fullfile(dataOutputFolder,sprintf('%s_%d_cim_output.mat',dataSource,ii));
        save(fOut,'R_cim');
        fOut = fullfile(dataOutputFolder,sprintf('%s_%d_knn1_output.mat',dataSource,ii));
        save(fOut,'R_knn1');
        fOut = fullfile(dataOutputFolder,sprintf('%s_%d_knn6_output.mat',dataSource,ii));
        save(fOut,'R_knn6');
        fOut = fullfile(dataOutputFolder,sprintf('%s_%d_knn20_output.mat',dataSource,ii));
        save(fOut,'R_knn20');
        fOut = fullfile(dataOutputFolder,sprintf('%s_%d_vme_output.mat',dataSource,ii));
        save(fOut,'R_vmeMI');
        fOut = fullfile(dataOutputFolder,sprintf('%s_%d_ap_output.mat',dataSource,ii));
        save(fOut,'R_apMI');
    end
end

