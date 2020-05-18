function [data, names] = loadSimulatedData(dirName,dataLines)

% CAUTION CAUTION CAUTION: TEMPORARY LOADER! BUGS ABOUND
dirFiles = dir(fullfile(dirName,'*.csv'));

data = zeros(length(dirFiles),6);
names = cell(length(dirFiles),1);


% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Sr", "k1", "k2", "beta", "sigma", "maxBOLD"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


for i = 1:length(dirFiles)
    data(i,:) = table2array(readtable(fullfile(dirFiles(i).folder,dirFiles(i).name),opts));
    names{i} = dirFiles(i).name;
end

numbers = [];
% sort data numerically by run (assuming random ('false') followed by qp
% ('true') in the simulations
for i = 1:length(names)/2
    temp = strsplit(names{i},{'_','.'});
    numbers(i) = str2double(temp(3));
end

numbers = numbers';
[~,randomIndex] = sort(numbers);

randomData = data(randomIndex,:);
randomNames = names(randomIndex);

qpIndex = randomIndex+50;
qpData = data(qpIndex,:);
qpNames = names(qpIndex);

% hardcoding data from test_1
realData(1:50,:) = repmat([1. .1 .2 1.00 .4 1.5],50,1);

mean(qpData)
mean(randomData)
realData(1,:)

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
realDataPlot = realData(first,:);
realDataPlot(4) = qpPlotParams(4);

realData(first,:);

figure
stimulusFreqHzFine = logspace(log10(.01),log10(64),100);
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,qpPlotParams),'-r');
hold on
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,randomPlotParams),'-b');
%semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,realDataPlot)-.036,'-k');
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,realDataPlot),'-g');
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,realData(first,:)),'-k');
legend('q+','random','Veridical, scaled & transposed','Veridical','Northwest');
