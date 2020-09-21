function [stats] = processSimulatedData(dirName,groupingVars)


resultsNames = dir(fullfile(dirName,'*.mat'));
% Check if data have already been processed.
if ~isempty(find(strcmp({resultsNames.name}, 'loadedData.mat')==1))
    load(fullfile(dirName,'loadedData.mat'));
else
    data = loadSimulatedData(dirName);
end

dataTable = struct2table(data.qpfmriResults);

for i = 1:height(dataTable)
    temp = dataTable.entropyOverTrials(i,end);
    dataTable.entropyFinal(i) = temp{end}(end);
    dataTable.maxBOLDFinalPreFit(i) = dataTable.maxBOLDoverTrials(i,end);
    if contains(dataTable.outNum{i},'qp')
        dataTable.qpPresFinal{i} = 'qp';
    else
        dataTable.qpPresFinal{i} = 'random';
    end
end



inputVars = {'nOutcomes','noiseSD','qpPresFinal','maxBOLDFinal',...,
    'maxBOLDFinalPreFit','psiParamsBadsFinal','entropyFinal'};

dsa = dataTable(:,inputVars);

stats = grpstats(dsa,groupingVars,{'mean','median','std','meanci'})