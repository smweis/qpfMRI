function [qpfmriResults] = compiledSimulate(modelName, varargin)
%% [qpfmriResults] = compiledSimulate(modelName, varargin)
% compiledSimulate formats CLI input for simulation simulate
% 
% 
% Description:
%	This is the compiled version of simulate. It is designed to work with a
%	high performance computing system, which can pass only strings as
%	arguments. See simulate.m for how optional arguments should be specified.
%
% Inputs:
%   modelName                 - String. Equivalent to a function handle.  
%   nParams                   - String. Number of parameters for the
%                               model.
%   NOTE: These variables should be passed for all parameters where X is
%   the number of parameters. Right now, we support up to 5 parameters. 
%
%   parameterXName            - String. Name of parameter X.
%   parameterXLower           - String. Lower bound of parameter X.
%   parameterXUpper           - String. Upper bound of parameter X.
%   paramXnDivisions          - String. Number of divisions between lower
%                               and upper bound for parameter.
%   parameterXSpacing         - String. Spacing of parameter X. Can be log,
%                               linear, or zeno's. (See makeDomain.m)
%
% See simulate.m for all optional inputs. 

%Example: 
%{
% Logisitic Regression example
modelName = 'logistic';
% How many parameters, and what are their names?
nP = '4';
p1Name = 'slope';
p2Name = 'semiSat';
p3Name = 'beta';
p4Name = 'sigma';

% Parameters lower bounds on domains
p1L = '-2'; 
p2L = '.01';
p3L = '.5';
p4L = '.3';

% Parameters upper bounds on domains
p1U = '0';
p2U = '1';
p3U = '1.5';
p4U = '2.0';

% Parameters nDivisions for domain creation
p1N = '10';
p2N = '10';
p3N = '11';
p4N = '8';

% Parameters spacing for domain creations
p1S = 'log';
p2S = 'lin';
p3S = 'zeno';
p4S = 'lin';

% Q+ control
qpPres = 'true';

% Stimulus domain
stimLower = '.01';
stimUpper = '1';
stimnDivisions = '30';
stimSpacing = 'lin';

% Show plots just for verification that this works the same with simulate. 
showPlots = true;

% Tiny noise so we know if things are working
noiseSD = '.01';

nTrials = '30';

p1Sim = '.12';
p2Sim = '.30';
p3Sim = '1.0';
p4Sim = '.1';

% Note, this will save a copy of questData after it is initialized. 
[qpfmriResults]=compiledSimulate(modelName,'nParams',nP,...,
'param1Name',p1Name,'param2Name',p2Name,'param3Name',p3Name,'param4Name',p4Name,...,
'param1Lower',p1L,'param1Upper',p1U,'param1nDivisions',p1N,'param1Spacing',p1S,...,
'param2Lower',p2L,'param2Upper',p2U,'param2nDivisions',p2N,'param2Spacing',p2S,...,
'param3Lower',p3L,'param3Upper',p3U,'param3nDivisions',p3N,'param3Spacing',p3S,...,
'param4Lower',p4L,'param4Upper',p4U,'param4nDivisions',p4N,'param4Spacing',p4S,...,
'stimDomainUpper',stimUpper,'stimDomainLower',stimLower,...,
'stimDomainnDivisions',stimnDivisions,'stimDomainSpacing',stimSpacing,...,
'showPlots',showPlots,'qpPres',qpPres,'noiseSD',noiseSD,'nTrials',nTrials,...,
'param1Simulated',p1Sim,'param2Simulated',p2Sim,'param3Simulated',p3Sim,'param4Simulated',p4Sim);


%}
%% TO DO
%1. We are very inflexible with the number of parameters we can handle at
%the moment, limited by the fact that we're expecting only strings as
%inputs. Could point the compiled code to a struct with necessary
%parameters?


%% Handle inputs
p = inputParser;

% Only MATLAB required input is the model name. But in fact, we also need
% all strings to create the paramDomains for simulate.m to run properly. 
p.addRequired(modelName,@ischar);

%% Required input for any model.
% The following are required to assemble params domain.
p.addParameter('nParams','',@ischar);
p.addParameter('param1Name','',@ischar);
p.addParameter('param2Name','',@ischar);
p.addParameter('param3Name','',@ischar);
p.addParameter('param4Name','',@ischar);
p.addParameter('param5Name','',@ischar);
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

% Spacing can be 'log', 'lin', or 'zeno' for the beta parameter. 
% Zeno will asymptotically approach the midpoint of upper and lower with
% (nDivisions-1)/2 for each side.
% Note, nDivisions for zeno-spaced parameters must be odd.
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

% The following are to assemble the stimulus domain.
p.addParameter('stimDomainLower','',@ischar);
p.addParameter('stimDomainnDivisions','',@ischar);
p.addParameter('stimDomainUpper','',@ischar);
p.addParameter('stimDomainSpacing','lin',@ischar);

% Optional params
p.addParameter('qpPres','false',@ischar);
p.addParameter('headroom', '0.1', @ischar);
p.addParameter('maxBOLDInitialGuess', '1.0', @ischar);
p.addParameter('maxBOLDSimulated', '1.5', @ischar);
p.addParameter('TR','800', @ischar);
p.addParameter('trialLength','12', @ischar);
p.addParameter('outNum','test',@ischar);
p.addParameter('outFolder','Results',@ischar);
p.addParameter('seed','choose');
p.addParameter('nTrials','10',@ischar);
p.addParameter('stimulusStructDeltaT','100',@ischar);
p.addParameter('baselineStimulus','');
p.addParameter('maxBOLDStimulus','');
p.addParameter('nOutcomes','15',@ischar);
p.addParameter('noiseSD','.1',@ischar);
p.addParameter('baselineMaxBOLDInitial','6',@ischar);
p.addParameter('baselineMaxBOLDRepeating','5',@ischar);

% Plotting parameters. This will only accept logicals (not strings) and
% this will not be able to be used in compiled form. For testing, though,
% you can view plots from simulate.
p.addParameter('showPlots',false,@islogical);

p.parse(modelName, varargin{:});

%% Create the model function from the model name
model = str2func(modelName);


%% Build parameter domains
% Initialize paramDomains
paramsDomain = struct;

% Assemble the parameter domains
for name = 1:str2double(p.Results.nParams)
    lower = str2double(p.Results.(['param' num2str(name) 'Lower']));
    nDivision = str2double(p.Results.(['param' num2str(name) 'nDivisions']));
    upper = str2double(p.Results.(['param' num2str(name) 'Upper']));
    % Call the function to make the parameterDomain space
    paramsDomain.(p.Results.(['param' num2str(name) 'Name'])) = makeDomain(lower,upper,nDivision,'spacing',p.Results.(['param' num2str(name) 'Spacing']));
end

%% Get the default qpfmri parameters (we'll add the settings in later)
[myQpfmriParams,myQpParams] = qpfmriParams(model,paramsDomain);

%% Build stimulus domain

% Assemble the stimulus domain
if ~isempty(p.Results.stimDomainLower)
    lower = str2double(p.Results.stimDomainLower);
    nDivision = str2double(p.Results.stimDomainnDivisions);
    upper = str2double(p.Results.stimDomainUpper);
    % Make the sure the parameter domain can be properly specified.
    myQpfmriParams.stimulusDomain = {makeDomain(lower,upper,nDivision,'spacing',p.Results.stimDomainSpacing)};
    fprintf('Stimulus domain added from input.\n');
end

% Check if a baseline stimulus was passed
if ~isempty(p.Results.baselineStimulus)
    myQpfmriParams.baselineStimulus = p.Results.baselineStimulus;
end
% Check if a maxBOLD stimulus was passed
if ~isempty(p.Results.maxBOLDStimulus)
    myQpfmriParams.maxBOLDStimulus = p.Results.maxBOLDStimulus;
end


%% Handle simulated psychomeric parameters

% Check if simulated parameters were passed. If not, simulate.m
% will choose one randomly from the parameter domains.
for name = 1:str2double(p.Results.nParams)
    if ~isempty(p.Results.(['param' num2str(name) 'Simulated']))
        param = str2double(p.Results.(['param' num2str(name) 'Simulated']));
        myQpfmriParams.simulatedPsiParams(name) = param;
        fprintf('Parameter %s specified as %0.5f\n',myQpfmriParams.paramNamesInOrder{name},param);
    else
        fprintf('Parameter %s randomly chosen as %0.5f',myQpfmriParams.paramNamesInOrder{name},myQpfmriParams.simulatedPsiParams(name));
    end
end

%% Handle the rest of the optional arguments: 
if strcmpi(p.Results.qpPres,'true')
    myQpfmriParams.qpPres = true;
else
    myQpfmriParams.qpPres = false;
end

myQpfmriParams.headroom = str2double(p.Results.headroom);
myQpfmriParams.maxBOLDInitialGuess = str2double(p.Results.maxBOLDInitialGuess);
myQpfmriParams.maxBOLDSimulated = str2double(p.Results.maxBOLDSimulated);
myQpfmriParams.TR = str2double(p.Results.TR);
myQpfmriParams.trialLength = str2double(p.Results.trialLength);
myQpfmriParams.outNum = p.Results.outNum;
myQpfmriParams.outFolder = p.Results.outFolder;
myQpfmriParams.nTrials = str2double(p.Results.nTrials);
myQpfmriParams.stimulusStructDeltaT = str2double(p.Results.stimulusStructDeltaT);
myQpfmriParams.nOutcomes = str2double(p.Results.nOutcomes);
myQpfmriParams.noiseSD = str2double(p.Results.noiseSD);
myQpfmriParams.fMRInoise = myQpfmriParams.noiseSD .* 2*myQpfmriParams.maxBOLDSimulated;
myQpfmriParams.baselineMaxBOLDInitial = str2double(p.Results.baselineMaxBOLDInitial);
myQpfmriParams.baselineMaxBOLDRepeating = str2double(p.Results.baselineMaxBOLDRepeating);

% Set seed, or let the simulation script choose. 
if strcmpi(p.Results.seed,'choose')
    myQpfmriParams.seed = p.Results.seed;
else
    myQpfmriParams.seed = str2double(p.Results.seed);
end

fprintf('seed = %d', myQpfmriParams.seed);

%% Run Simulation
[qpfmriResults]=simulate(myQpfmriParams, myQpParams,'showPlots',p.Results.showPlots);

end


