function [data] = loadSimulatedData(dirName,dataLines)

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
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
optsResults.VariableNames = ["Sr", "k1", "k2", "beta", "sigma", "maxBOLD", "qpPres","simID"];
optsResults.VariableTypes = ["double", "double", "double", "double", "double", "double", "string", "string"];

optsParams.VariableNames = ["SrSim","k1Sim","k2Sim","betaSim","sigmaSim","maxBOLDSim","simID"];
optsParams.VariableTypes = ["double", "double", "double", "double", "double", "double", "string"];

paramsData = table('Size',[length(resultsDirFiles),length(optsParams.VariableNames)],...
    'VariableTypes',optsParams.VariableTypes,'VariableNames',optsParams.VariableNames);

resultsData = table('Size',[length(resultsDirFiles),length(optsResults.VariableNames)],...
    'VariableTypes',optsResults.VariableTypes,'VariableNames',optsResults.VariableNames);

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


for i = 1:length(resultsDirFiles)
    temp = strsplit(resultsDirFiles(i).name,{'_','.'});
    paramsData.simID(i) = [temp{1} '_' temp{2}];
    resultsData.qpPres(i) = temp(1);
    resultsData.simID(i) = [temp{1} '_' temp{2}];
    resultsData(i,1:6) = readtable(fullfile(resultsDirFiles(i).folder,resultsDirFiles(i).name));
    paramsData(i,1:6) = readtable(fullfile(paramDirFiles(i).folder,paramDirFiles(i).name));
end

data = join(resultsData,paramsData);

avgResults = varfun(@mean,data,'InputVariables',{'Sr','k1','k2','beta','sigma'},...
       'GroupingVariables','qpPres')

   

% If params are all the same, we can use this: 
sampleSimulatedParams = [data.SrSim(1) data.k1Sim(1) data.k2Sim(1) data.betaSim(1),data.sigmaSim(1)];

idx = avgResults.qpPres=="true";
a = avgResults(idx,:);
qpParams = [a.mean_Sr a.mean_k1 a.mean_k2 a.mean_beta a.mean_sigma];
idx = avgResults.qpPres=="false";
b = avgResults(idx,:);
randomParams = [b.mean_Sr b.mean_k1 b.mean_k2 b.mean_beta b.mean_sigma];


% Plot the average results
figure
freqDomain = logspace(log10(.01),log10(100),100);
predictedRelativeResponse = doeTemporalModel(freqDomain,sampleSimulatedParams) - ...
        doeTemporalModel(0,sampleSimulatedParams);
predictedQpRelativeResponse = doeTemporalModel(freqDomain,qpParams) - ...
        doeTemporalModel(0,qpParams);
predictedRandomRelativeResponse = doeTemporalModel(freqDomain,randomParams) - ...
        doeTemporalModel(0,randomParams);

semilogx(freqDomain,predictedRelativeResponse,'-k','LineWidth',2);
hold on
semilogx(freqDomain,predictedQpRelativeResponse,'-r','LineWidth',2);
semilogx(freqDomain,predictedRandomRelativeResponse,'-b','LineWidth',2);

set(gca,'XScale', 'log')
legend('Veridical','Q+','Random','Location','Northwest');



data.SrBias = (data.Sr - data.SrSim)./data.SrSim;
data.k1Bias = (data.k1 - data.k1Sim)./data.k1Sim;
data.k2Bias = (data.k2 - data.k2Sim)./data.k2Sim;
data.betaBias = (data.beta - data.betaSim)./data.betaSim;
data.sigmaBias = (data.sigma - data.sigmaSim)./data.sigmaSim;


data.SrNoise = abs(data.Sr - data.SrSim)./data.SrSim;
data.k1Noise = abs(data.k1 - data.k1Sim)./data.k1Sim;
data.k2Noise = abs(data.k2 - data.k2Sim)./data.k2Sim;
data.betaNoise = abs(data.beta - data.betaSim)./data.betaSim;
data.sigmaNoise = abs(data.sigma - data.sigmaSim)./data.sigmaSim;


biasResults = varfun(@mean,data,'InputVariables',{'SrBias','k1Bias','k2Bias','betaBias','sigmaBias'},...
       'GroupingVariables','qpPres')

   
noiseResults = varfun(@mean,data,'InputVariables',{'SrNoise','k1Noise','k2Noise','betaNoise','sigmaNoise'},...
       'GroupingVariables','qpPres')

