%% Script which processes the netbenchmark data
clear;
clc;

% setup and start the parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 6;
saveProfile(myCluster);
p = gcp

if(ispc)
    dataRepo = 'C:\\Users\\Kiran\\Documents\\data\\netbenchmark\\inputs';
    dataOutputFolder = 'C:\\Users\\Kiran\\Documents\\data\\netbenchmark\\matlab_outputs';
else
    dataRepo = '/data/netbenchmark/inputs/';
    dataOutputFolder = '/data/netbenchmark/matlab_outputs';
end

dataSourcesToProcess = {'syntren300'};
numDatasets = 200;

globalNoiseVec = [0,10,20,30,40,50];
localNoiseVec  = [0,10,20,30,40,50];

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for dataSourceIdx=1:length(dataSourcesToProcess)
    dataSource = dataSourcesToProcess{dataSourceIdx};
    for gn=globalNoiseVec
        for ln=localNoiseVec
            for ii=1:numDatasets
                dispstat(sprintf('Computing %s{%d,%d}[%d/%d]',dataSource, gn, ln, ii, numDatasets),'keepthis', 'timestamp');
                subfolder = sprintf('gn_%d_ln_%d',gn,ln);
                fIn = fullfile(dataRepo,subfolder,sprintf('%s_%d.mat',dataSource,ii));
                load(fIn);
                
                cimfOut = fullfile(dataOutputFolder,subfolder,sprintf('%s_%d_cim_output.mat',dataSource,ii));
                knn1fOut = fullfile(dataOutputFolder,subfolder,sprintf('%s_%d_knn1_output.mat',dataSource,ii));
                knn6fOut = fullfile(dataOutputFolder,subfolder,sprintf('%s_%d_knn6_output.mat',dataSource,ii));
                knn20fOut = fullfile(dataOutputFolder,subfolder,sprintf('%s_%d_knn20_output.mat',dataSource,ii));
                vmefOut = fullfile(dataOutputFolder,subfolder,sprintf('%s_%d_vme_output.mat',dataSource,ii));
                apfOut = fullfile(dataOutputFolder,subfolder,sprintf('%s_%d_ap_output.mat',dataSource,ii));
                taufOut = fullfile(dataOutputFolder,subfolder,sprintf('%s_%d_tau_output.mat',dataSource,ii));
                
                if(~exist(cimfOut,'file') || ~exist(knn1fOut,'file') || ~exist(knn6fOut,'file') || ...
                   ~exist(knn20fOut,'file') || ~exist(vmefOut,'file') || ~exist(apfOut,'file') || ~exist(taufOut,'file'))
                    X = data;
                    [R_cim,R_knn1,R_knn6,R_knn20,R_vmeMI,R_apMI,R_tau] = pair_all_cc_mex( X );

                    if ~exist(fullfile(dataOutputFolder,subfolder), 'dir')
                        mkdir(fullfile(dataOutputFolder,subfolder));
                    end

                    save(cimfOut,'R_cim');
                    save(knn1fOut,'R_knn1');
                    save(knn6fOut,'R_knn6');
                    save(knn20fOut,'R_knn20');
                    save(vmefOut,'R_vmeMI');
                    save(apfOut,'R_apMI');
                    save(taufOut,'R_tau');
                end
            end
        end
    end
end