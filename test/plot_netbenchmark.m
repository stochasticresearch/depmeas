%% Plot the netbenchmark combined data
clear;
clc;
close all;

datasetsToPlot = {'syntren300'};
globalNoiseVec = [0,10,20,30,40,50];
localNoiseVec = [10,20,30,40,50];
if(ispc)
    combinedDataRepo = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\netbenchmark';
elseif(ismac)
    combinedDataRepo = '/data/netbenchmark/combined_outputs';
else
    combinedDataRepo = '/data/netbenchmark/combined_outputs';
end

%estimatorsToPlot = {'"CIM"', '"Kendall"', '"KNN1"', '"KNN6"', '"KNN20"', '"vME"', '"AP"', '"rand"'};
estimatorsToPlot = {'"CIM"', '"KNN1"', '"KNN6"', '"KNN20"', '"vME"', '"AP"'};
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
fontSizeVal = 20;
for datasetIdx=1:length(datasetsToPlot)
    for gn=globalNoiseVec
        for ln=localNoiseVec
            % if the file exists -- process it
            ds = datasetsToPlot{datasetIdx};
            combinedDataFile = fullfile(combinedDataRepo,sprintf('%s_combined_gn_%d_ln_%d.csv',ds,gn,ln));
            if(exist(combinedDataFile,'file'))
                figure;
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
                title(sprintf('%s - [%d/%d]',ds,gn,ln), 'FontSize', fontSizeVal);
                grid on;
                ylabel('max[AUPR_{20}]', 'FontSize', fontSizeVal);
                hh = gca;
                hh.FontSize = fontSizeVal;
                set(h,'LineWidth',3);
                set(gca, 'YScale', 'log');
            end
        end
    end
end