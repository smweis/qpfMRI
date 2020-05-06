function [data, names] = loadSimulatedData(dirName,dataLines)

% CAUTION CAUTION CAUTION: TEMPORARY LOADER! BUGS ABOUND

dirFiles = dir(fullfile('..',dirName,'*.csv'));

data = zeros(length(dirFiles),5);
names = cell(length(dirFiles),1);


% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"], "ThousandsSeparator", ",");


for i = 1:length(dirFiles)
    data(i,:) = table2array(readtable(fullfile(dirFiles(i).folder,dirFiles(i).name),opts));
    names{i} = dirFiles(i).name;
end


% TO DO: SORT output numerically. 
%1. 

