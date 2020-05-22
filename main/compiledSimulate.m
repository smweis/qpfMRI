function [psiParamsFit,maxBOLD,questDataCopy] = compiledSimulate(modelName, varargin)
%% compiledSimulate formats CLI input for simulation
%% simulate
% A script that will simulate fMRI BOLD data and fit a model with or
% without Q+ control. 
%
% 
% Description:
%	Takes in a model and a possible set of parameters and whether or not Q+ 
%   is in control of things flag. 
%
% Inputs:
%   model                 - String that is equivalent to a function handle. 
%                           This should be the 'continuous function.' 
%
%
% See simulate.m for all additional inputs. 

%Example: 
%{

modelName = 'doeTemporalModel';

% Parameters 1-5 lower bounds on domains
p1L = '.899'; 
p2L = '.01';
p3L = '.01';
p4L = '.8';
p5L = '.3';
% Parameters 1-5 upper bounds on domains
p1U = '1.099';
p2U = '.4';
p3U = '.4';
p4U = '1.4';
p5U = '1.0';
% Parameters 1-5 intervals for domain creation
p1I = '.025';
p2I = '.04';
p3I = '.04';
p4I = '.1';
p5I = '.2';

% Simulated values for parameters. (These are optional).
p1S = '1.004';
p2S = '.016';
p3S = '.118';
p4S = '1.0';
p5S = '.1';

% Show plots just for verification that this works the same with simulate. 
showPlots = true;

% Note, this will save a copy of questData after it is initialized. 
[psiParamsFit,maxBOLD,questDataCopy]=compiledSimulate(modelName,...,
'param1Lower',p1L,'param1Upper',p1U,'param1nDivisions',p1I,...,
'param2Lower',p2L,'param2Upper',p2U,'param2nDivisions',p2I,...,
'param3Lower',p3L,'param3Upper',p3U,'param3nDivisions',p3I,...,
'param4Lower',p4L,'param4Upper',p4U,'param4nDivisions',p4I,...,
'param5Lower',p5L,'param5Upper',p5U,'param5nDivisions',p5I,...,
'param1Simulated',p1S,'param2Simulated',p2S,'param3Simulated',p3S,...,
'param4Simulated',p4S,'param5Simulated',p5S,'showPlots',showPlots);


% Logisitic Regression example
modelName = 'logistic';

% Parameters lower bounds on domains
p1L = '.01'; 
p2L = '.01';
p3L = '.8';
p4L = '.3';

% Parameters upper bounds on domains
p1U = '1';
p2U = '1';
p3U = '1.4';
p4U = '2.0';

% Parameters nDivisions for domain creation
p1N = '30';
p2N = '30';
p3N = '8';
p4N = '8';

% Parameters spacing for domain creations
p1S = 'lin';
p2S = 'lin';
p3S = 'lin';
p4S = 'lin';

% NoiseSD sampling
noiseSDLower = '.1';
noiseSDUpper = '.8';
noiseSDnDivisions = '7';

% Q+ control
qpPres = 'true';

% Stimulus domain
stimLower = '.01';
stimUpper = '1';
stimnDivisions = '30';
stimSpacing = 'lin';

% Show plots just for verification that this works the same with simulate. 
showPlots = true;

% Note, this will save a copy of questData after it is initialized. 
[psiParamsFit,maxBOLD,questDataCopy]=compiledSimulate(modelName,...,
'param1Lower',p1L,'param1Upper',p1U,'param1nDivisions',p1N,'param1Spacing',p1S,...,
'param2Lower',p2L,'param2Upper',p2U,'param2nDivisions',p2N,'param2Spacing',p2S,...,
'param3Lower',p3L,'param3Upper',p3U,'param3nDivisions',p3N,'param3Spacing',p3S,...,
'param4Lower',p4L,'param4Upper',p4U,'param4nDivisions',p4N,'param4Spacing',p4S,...,
'noiseSDLower',noiseSDLower,'noiseSDUpper',noiseSDUpper,'noiseSDnDivisions',noiseSDnDivisions,...,
'stimDomainUpper',stimUpper,'stimDomainLower',stimLower,...,
'stimDomainnDivisions',stimnDivisions,'stimDomainSpacing',stimSpacing,...,
'showPlots',showPlots,'qpPres',qpPres);


%}


%% Handle inputs
p = inputParser;

% Only MATLAB required input is the model name. But in fact, we also need
% all strings to create the paramDomains for simulate.m to run properly. 
p.addRequired(modelName,@ischar);

%% Required input for any model.
% The following are required to assemble params domain.
p.addParameter('param1Lower','',@ischar);
p.addParameter('param2Lower','',@ischar);
p.addParameter('param3Lower','',@ischar);
p.addParameter('param4Lower','',@ischar);
p.addParameter('param5Lower','',@ischar);
p.addParameter('param1nDivisions','',@ischar);
p.addParameter('param2nDivisions','',@ischar);
p.addParameter('param3nDivisions','',@ischar);
p.addParameter('param4nDivisions','',@ischar);
p.addParameter('param5nDivisions','',@ischar);
p.addParameter('param1Upper','',@ischar);
p.addParameter('param2Upper','',@ischar);
p.addParameter('param3Upper','',@ischar);
p.addParameter('param4Upper','',@ischar);
p.addParameter('param5Upper','',@ischar);
p.addParameter('param1Spacing','lin',@ischar);
p.addParameter('param2Spacing','lin',@ischar);
p.addParameter('param3Spacing','lin',@ischar);
p.addParameter('param4Spacing','lin',@ischar);
p.addParameter('param5Spacing','lin',@ischar);

% The following are optional for specific simulated parameters. 
% Otherwise parameters will be selected randomly from the parameter domain.
p.addParameter('param1Simulated','',@ischar);
p.addParameter('param2Simulated','',@ischar);
p.addParameter('param3Simulated','',@ischar);
p.addParameter('param4Simulated','',@ischar);
p.addParameter('param5Simulated','',@ischar);

% To provide a range for noiseSD, enter the same values:
p.addParameter('noiseSDLower','.1',@ischar);
p.addParameter('noiseSDUpper','.1',@ischar);
p.addParameter('noiseSDnDivisions','1',@ischar);

% The following are to assemble the stimulus domain.
p.addParameter('stimDomainLower','',@ischar);
p.addParameter('stimDomainnDivisions','',@ischar);
p.addParameter('stimDomainUpper','',@ischar);
p.addParameter('stimDomainSpacing','lin',@ischar);

%% Optional model general parameters
p.addParameter('qpPres','false',@ischar);
p.addParameter('nTrials','30',@ischar);
p.addParameter('stimulusStructDeltaT','100',@ischar);
p.addParameter('maxBOLDSimulated','1.6',@ischar);
p.addParameter('maxBOLD','1.0',@ischar);
p.addParameter('baselineStimulus','',@ischar);
p.addParameter('maxBOLDStimulus','',@ischar);
p.addParameter('nOutcomes','51',@ischar);
p.addParameter('headroom','.1',@ischar);
p.addParameter('TR','800', @ischar);
p.addParameter('trialLength','12', @ischar);
p.addParameter('outNum','test',@ischar);
p.addParameter('seed','choose',@ischar);

% Optional params for future support (not currently supported)
p.addParameter('questDataCopy',{},@isstruct);

% Plotting parameters. This will only accept logicals (not strings) and
% this will not be able to be used in compiled form. For testing, though,
% you can view plots from simulate.
p.addParameter('showPlots',false,@islogical);

p.parse(modelName, varargin{:});

%% Create the model function from the model name
model = str2func(modelName);
[paramNamesInOrder] = checkModel(model);

%% Build parameter domains
% Assemble the paramsDomain and, if specified, simulatedPsiParams

% Initialize paramDomains
paramsDomain = struct;

% Assemble the parameter domains
for name = 1:length(paramNamesInOrder)
    lower = str2double(p.Results.(['param' num2str(name) 'Lower']));
    nDivision = str2double(p.Results.(['param' num2str(name) 'nDivisions']));
    upper = str2double(p.Results.(['param' num2str(name) 'Upper']));
    % Make the sure the parameter domain can be properly specified.
    assert(lower<upper,'Improper parameter domain %s',paramNamesInOrder{name});
    if contains(p.Results.(['param' num2str(name) 'Spacing']),'lin')
        paramsDomain.(paramNamesInOrder{name}) = linspace(lower,upper,nDivision);
    elseif contains(p.Results.(['param' num2str(name) 'Spacing']),'log')
        paramsDomain.(paramNamesInOrder{name}) = logspace(lower,upper,nDivision);
    else
        error('Improper spacing key for %s. Log or linear spacing supported',paramNamesInOrder{name});
    end
end

%% Build stimulus domain
% Initialize 
stimulusDomain = {};

% Assemble the parameter domains
if ~isempty(p.Results.stimDomainLower)
    lower = str2double(p.Results.stimDomainLower);
    nDivision = str2double(p.Results.stimDomainnDivisions);
    upper = str2double(p.Results.stimDomainUpper);
    % Make the sure the parameter domain can be properly specified.
    assert(lower<upper,'Improper parameter domain %s',paramNamesInOrder{name});
    if contains(p.Results.stimDomainSpacing,'lin')
        stimulusDomain = {linspace(lower,upper,nDivision)};
    elseif contains(p.Results.stimDomainSpacing,'log')
        stimulusDomain = {logspace(lower,upper,nDivision)};
    else
        error('Improper spacing key for stimulus domain. Log or linear spacing supported');
    end

end




% Check if one of the simulated parameters was passed. If not, simulate.m
% will choose one randomly from the parameter domains.
if ~isempty(p.Results.param1Simulated)
    simulatedPsiParams = struct;
    for name = 1:length(paramNamesInOrder)
        param = str2double(p.Results.(['param' num2str(name) 'Simulated']));
        simulatedPsiParams.(paramNamesInOrder{name}) = param;
    end
else
    simulatedPsiParams = {};
end

% To provide a range for noiseSD, enter the same values:
noiseSD = linspace(str2double(p.Results.noiseSDLower),...,
    str2double(p.Results.noiseSDUpper),str2double(p.Results.noiseSDnDivisions));


% Handle the rest of the optional arguments: 
if strcmpi(p.Results.qpPres,'true')
    qpPres = true;
else
    qpPres = false;
end

nTrials = str2double(p.Results.nTrials);
stimulusStructDeltaT = str2double(p.Results.stimulusStructDeltaT);
maxBOLDSimulated = str2double(p.Results.maxBOLDSimulated);
maxBOLD = str2double(p.Results.maxBOLD);
baselineStimulus = p.Results.baselineStimulus;
maxBOLDStimulus = p.Results.maxBOLDStimulus;
nOutcomes = str2double(p.Results.nOutcomes);
headroom = str2double(p.Results.headroom);
TR = str2double(p.Results.TR);
trialLength = str2double(p.Results.trialLength);
outNum = p.Results.outNum;
seed = str2double(p.Results.seed);

[psiParamsFit,maxBOLD,questDataCopy] = simulate(model, paramsDomain,...,
    'simulatedPsiParams',simulatedPsiParams,'qpPres',qpPres,...,
    'nTrials',nTrials,'stimulusStructDeltaT',stimulusStructDeltaT,...,
    'maxBOLDSimulated',maxBOLDSimulated,'maxBOLD',maxBOLD,...,
    'baselineStimulus',baselineStimulus,'maxBOLDStimulus',maxBOLDStimulus,...,
    'nOutcomes',nOutcomes,'headroom',headroom,'noiseSD',noiseSD,...,
    'TR',TR,'trialLength',trialLength,'outNum',outNum,'seed',seed,...,
    'showPlots',p.Results.showPlots,'stimulusDomain',stimulusDomain);

end


