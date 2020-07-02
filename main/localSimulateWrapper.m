
% Provide a model handle
model = @logistic;

% Specify the parameter domain. Each value must correspond to a parameter
% expected by the model. 
paramsDomain = struct;
paramsDomain.slope = makeDomain(-1.2,-.2,10,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.beta = makeDomain(.75,1.25,11,'spacing','zeno');

% Sigma in the parameter domain is searching for noiseSD
paramsDomain.sigma = makeDomain(.5,4,8);

% Specify a stimulus domain and whether it spaced linear or log.
stimulusDomain = {makeDomain(.01,1,25)};
stimulusDomainSpacing = 'lin';

% Number of trials to run.
nTrials = 30;
trialLength = 12;
% Allow Q+ to control the stimuli or not (false).
qpPres = true;

nOutcomes = 5;% Set the number of outcome categories / bins.

showPlots = false; % Do you want to see plots?

simulatedPsiParams = struct;
simulatedPsiParams.slope = .35;
simulatedPsiParams.semiSat = .63;
simulatedPsiParams.beta = 1.0;

maxBOLDSimRange = makeDomain(.5,2,10);
noiseSDRange = makeDomain(.05,.3,10);
for i = 21:40
    
    maxBOLDSimulated = randsample(maxBOLDSimRange,1);
    noiseSD = randsample(noiseSDRange,1);
    outNum = strcat('qpControl_',num2str(i),'_noise_',num2str(noiseSD),'_');
    
    [psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain,...,
    'qpPres',qpPres, 'showPlots',showPlots,'stimulusDomain',stimulusDomain,...,
    'stimulusDomainSpacing',stimulusDomainSpacing,'noiseSD',noiseSD,...,
    'simulatedPsiParams',simulatedPsiParams,'nTrials',nTrials,...,
    'maxBOLDSimulated',maxBOLDSimulated,'trialLength',trialLength,'nOutcomes',nOutcomes,...,
    'outNum',outNum);
end
sameParams = true;
baseline = .01;
dirName = ('');
factorName = 'qpPres';
[data] = summarizeAndPlotSimulations(model,paramsDomain,factorName,...,
    sameParams,{stimulusDomain},baseline,dirName);