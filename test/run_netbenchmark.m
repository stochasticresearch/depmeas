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
numDatasets = 150;

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

                X = data;
                [R_cim,R_knn1,R_knn6,R_knn20,R_vmeMI,R_apMI,R_tau] = pair_all_cc_mex( X );

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
                fOut = fullfile(dataOutputFolder,sprintf('%s_%d_tau_output.mat',dataSource,ii));
                save(fOut,'R_tau');
            end
        end
    end
end