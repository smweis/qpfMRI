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


% sort data numerically by run (assuming random ('false') followed by qp
% ('true') in the simulations
for i = 1:length(names)/2
    temp = strsplit(names{i},{'_','.'});
    numbers(i) = str2double(temp(3));
end

numbers = numbers';
[~,randomIndex] = sort(numbers);

randomData = data(randomIndex,:);
randomNames = names(idx);

qpIndex = randomIndex+300;
qpData = data(qpIndex,:);
qpNames = names(qpIndex);

% hardcoding data from test_1
realData(1:50,:) = repmat([.98 .003 .06 1.00 .4],50,1);
realData(51:100,:) = repmat([.98 .003 .06 1.00 .8],50,1);
realData(101:150,:) = repmat([.92 .003 .06 1.00 .4],50,1);
realData(151:200,:) = repmat([.92 .007 .06 1.00 .4],50,1);
realData(201:250,:) = repmat([.92 .003 .12 1.00 .4],50,1);
realData(251:300,:) = repmat([.92 .003 .18 1.00 .4],50,1);


qpSigned = (qpData - realData)./realData;
randomSigned = (randomData - realData)./realData;

mean(qpSigned)
mean(randomSigned)

qpUnsigned = abs(qpData - realData)./realData;
randomUnsigned = abs(randomData - realData)./realData;

mean(qpUnsigned)
mean(randomUnsigned)


% one run and plotting
first = 1;
last = first+49;

qpOneRunData = qpData(first:last,:);
randomOneRunData = randomData(first:last,:);

qpPlotParams = mean(qpOneRunData);
randomPlotParams = mean(randomOneRunData);

% Hold beta at 1 for plotting
qpPlotParams(4) = 1;
randomPlotParams(4) = 1;

realData(first,:);

figure
stimulusFreqHzFine = logspace(log10(.01),log10(64),100);
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,qpPlotParams),'-r');
hold on
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,randomPlotParams),'-b');
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,realData(first,:)),'-k');
