function [subject,myQpfmriParams,myQpParams] = getparams()

% Load global parameters
fid = fopen("C:\Users\jacob.frank\Documents\blue\rtQuest\derivatives\realTime\test\processed\run0\global.json");
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
global_params = jsondecode(str);

subject = global_params.subject;
run = global_params.run;

% Provide a model handle
model = @logistic;

% Specify the parameter domain. Each value must correspond to a parameter
% expected by the model. 

paramsDomain = struct;
paramsDomain.slope = makeDomain(-1.2,-.2,10,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.beta = makeDomain(.75,1.25,11,'spacing','zeno');

% Sigma in the parameter domain is searching for noiseSD
paramsDomain.sigma = makeDomain(0.01,1,10);

% Specify a stimulus domain and whether it spaced linear or log.
stimulusDomain = {makeDomain(.5,1,12)};
stimulusDomainSpacing = 'log';

% Number of trials to run and TR.
nTrials = 20;
TRmsecs = 800;

% Allow Q+ to control the stimuli or not (false).
qpPres = true;

% Set the number of outcome categories / bins.
nOutcomes = 15;

% Do you want to see plots?
%showPlots = true; 

%How long the trials are (in seconds).
trialLength = 12;

outNum = str2double(run);

[myQpfmriParams,myQpParams] = qpfmriParams(model,paramsDomain,'qpPres',qpPres,...,
'stimulusDomain',stimulusDomain,'stimulusDomainSpacing',stimulusDomainSpacing,...,
'nTrials',nTrials,'trialLength',trialLength,'nOutcomes',nOutcomes,'TR',TRmsecs,'outNum',outNum);

end


