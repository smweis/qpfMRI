function [data,myQpfmriParams] = loadSimulatedData(dirName)

% Grab all results .mat files from the directory dirName
resultsNames = dir(fullfile(dirName,'*.mat'));

% Check if data have already been processed.
if ~isempty(find(strcmp({resultsNames.name}, 'loadedData.mat')==1,1))
    warning('Data already loaded. Will reload data.');
    resultsNames = resultsNames(~startsWith({resultsNames.name}, 'loadedData'));
end

% Load first struct in.
data = load(fullfile(dirName,resultsNames(1).name));

% Combine all other structs into one struct. 
for i = 2:length(resultsNames)
    temp = load(fullfile(dirName,resultsNames(i).name));
    data = cell2struct(cellfun(@vertcat,struct2cell(data),struct2cell(temp),'uni',0),fieldnames(data),1);
end



% Save the file to the directory for easier loading.
save(fullfile(dirName,'loadedData'),'data');



