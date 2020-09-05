function [data] = loadSimulatedData(dirName,model, dataLines)

paramNamesInOrder = checkModel(model);

% CAUTION CAUTION CAUTION: TEMPORARY LOADER! BUGS ABOUND
paramDir = [dirName filesep 'params'];
resultsDir = [dirName filesep 'results'];


resultsDirFiles = dir(fullfile(resultsDir,'*.csv'));
paramDirFiles = dir(fullfile(paramDir,'*.csv'));


% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", length(paramNamesInOrder)+1);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
optsResults.VariableNames = ["slope", "semiSat", "beta", "sigma", "maxBOLD","simID"];
optsResults.VariableTypes = ["double", "double", "double", "double", "double", "string"];

optsParams.VariableNames = ["slopeSim","semiSatSim","betaSim","sigmaSim","maxBOLDSim","nOutcomes","simIDSim"];
optsParams.VariableTypes = ["double", "double", "double", "double", "double","double","string"];

paramsData = table('Size',[length(resultsDirFiles),length(optsParams.VariableNames)],...
    'VariableTypes',optsParams.VariableTypes,'VariableNames',optsParams.VariableNames);

resultsData = table('Size',[length(resultsDirFiles),length(optsResults.VariableNames)],...
    'VariableTypes',optsResults.VariableTypes,'VariableNames',optsResults.VariableNames);

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";



for i = 1:length(resultsDirFiles)
    temp = strsplit(resultsDirFiles(i).name,{'_','.'});
    paramsData.simIDSim(i) = [temp{1} '_' temp{2}];
    resultsData.qpPres(i) = temp(1);
    resultsData.simID(i) = [temp{1} '_' temp{2}];
    
    if contains(resultsDirFiles(i).name,'nOutcomes')
        resultsData.nOutcomes(i) = str2double(temp{4});
    end
    
    resultsData(i,1:length(paramNamesInOrder)+1) = readtable(fullfile(resultsDirFiles(i).folder,resultsDirFiles(i).name));
    paramsData(i,1:length(paramNamesInOrder)+2) = readtable(fullfile(paramDirFiles(i).folder,paramDirFiles(i).name));
end

sortrows(resultsData,'simID');
sortrows(paramsData,'simIDSim');
data = [resultsData paramsData];


%% Check the loaded data for duplicate values. 
% We want to check for duplicate values in our fit parameters. These should
% almost NEVER result in identical parameters; otherwise, we suspect
% something went wrong with setting the random seed. 
[~, dupeInd] = unique(data(:,1:length(paramNamesInOrder)+1), 'rows');
duplicate_ind = setdiff(1:size(data, 1:length(paramNamesInOrder)+1), dupeInd);
if ~isempty(duplicate_ind)
    warning('Duplicate values detected. Inspect data and logs to ensure unique simulations.');
end

 
