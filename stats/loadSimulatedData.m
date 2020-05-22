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

optsParams.VariableNames = ["slopeSim","semiSatSim","betaSim","sigmaSim","maxBOLDSim","simIDSim"];
optsParams.VariableTypes = ["double", "double", "double", "double", "double", "string"];

paramsData = table('Size',[length(resultsDirFiles),length(optsParams.VariableNames)],...
    'VariableTypes',optsParams.VariableTypes,'VariableNames',optsParams.VariableNames);

resultsData = table('Size',[length(resultsDirFiles),length(optsResults.VariableNames)],...
    'VariableTypes',optsResults.VariableTypes,'VariableNames',optsResults.VariableNames);

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


for i = 1:length(resultsDirFiles)
    temp = strsplit(resultsDirFiles(i).name,{'_','.'});
    paramsData.simIDSim(i) = [temp{1} '_' temp{3}];
    resultsData.qpPres(i) = temp(1);
    resultsData.simID(i) = [temp{1} '_' temp{3}];
    resultsData(i,1:length(paramNamesInOrder)+1) = readtable(fullfile(resultsDirFiles(i).folder,resultsDirFiles(i).name));
    paramsData(i,1:length(paramNamesInOrder)+1) = readtable(fullfile(paramDirFiles(i).folder,paramDirFiles(i).name));
end

sortrows(resultsData,'simID');
sortrows(paramsData,'simIDSim');
data = [resultsData paramsData];

avgResults = varfun(@mean,data,'InputVariables',paramNamesInOrder,...
       'GroupingVariables','qpPres')

   

% If params are all the same, we can use this: 
sampleSimulatedParams = [data.slopeSim(1) data.semiSatSim(1) data.betaSim(1),data.sigmaSim(1)]

idx = avgResults.qpPres=="true";
a = avgResults(idx,:);
qpParams = [a.mean_slope a.mean_semiSat a.mean_beta a.mean_sigma];
idx = avgResults.qpPres=="false";
b = avgResults(idx,:);
randomParams = [b.mean_slope b.mean_semiSat b.mean_beta b.mean_sigma];


% Plot the average results
figure
stimDomain = linspace(.01,1,20);
predictedRelativeResponse = model(stimDomain,sampleSimulatedParams) - ...
        model(0.01,sampleSimulatedParams);
predictedQpRelativeResponse = model(stimDomain,qpParams) - ...
        model(0.01,qpParams);
predictedRandomRelativeResponse = model(stimDomain,randomParams) - ...
        model(0,randomParams);

plot(stimDomain,predictedRelativeResponse,'-k','LineWidth',2);
hold on
plot(stimDomain,predictedQpRelativeResponse,'-r','LineWidth',2);
plot(stimDomain,predictedRandomRelativeResponse,'-b','LineWidth',2);

set(gca,'XScale', 'lin')
legend('Veridical','Q+','Random','Location','Northwest');
title(horzcat('Average across all from: ', dirName));


data.slopeBias = (data.slope - data.slopeSim)./data.slopeSim;
data.semiSatBias = (data.semiSat - data.semiSatSim)./data.semiSatSim;
data.betaBias = (data.beta - data.betaSim)./data.betaSim;
data.sigmaBias = (data.sigma - data.sigmaSim)./data.sigmaSim;


data.slopeNoise = abs(data.slope - data.slopeSim)./data.slopeSim;
data.semiSatNoise = abs(data.semiSat - data.semiSatSim)./data.semiSatSim;
data.betaNoise = abs(data.beta - data.betaSim)./data.betaSim;
data.sigmaNoise = abs(data.sigma - data.sigmaSim)./data.sigmaSim;


biasResults = varfun(@mean,data,'InputVariables',{'slopeBias','semiSatBias','betaBias','sigmaBias'},...
       'GroupingVariables','qpPres')

   
noiseResults = varfun(@mean,data,'InputVariables',{'slopeNoise','semiSatNoise','betaNoise','sigmaNoise'},...
       'GroupingVariables','qpPres')

   
function plotOneRun(data,row)

figure;
stimDomain = linspace(.01,1,100);

qpParams = [data.slope(row) data.semiSat(row) data.beta(row)];
sampleSimulatedParams = [data.slopeSim(row) data.semiSatSim(row) data.betaSim(row)];

predictedRelativeResponse = logistic(stimDomain,sampleSimulatedParams) - ...
        logistic(0.01,sampleSimulatedParams);
predictedQpRelativeResponse = logistic(stimDomain,qpParams) - ...
        logistic(0.01,qpParams);

plot(stimDomain,predictedRelativeResponse,'-k','LineWidth',2);
hold on;
plot(stimDomain,predictedQpRelativeResponse,'-r','LineWidth',2);

if contains(data.simID(row),'true')
    qpPres = 'Q+ stimulus selection';
else
    qpPres = 'Random stimulus selection';
end

legend('Veridical',qpPres,'Location','Northwest');

titleTxt = strcat('Single simulation: ',qpPres,' data row: ', num2str(row));
title(titleTxt);
hold off;

end


plotOneRun(data,1);

end   
