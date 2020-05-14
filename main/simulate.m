function [psiParamsFit,maxBOLD]=simulate(model, paramsDomain, qpPres, varargin)

% TO DO: HUNT FOR ANYWHERE IN THE SCRIPT THAT REFERENCES A PARAMETER BY
% NUMBER RATHER THAN FINDING IT BY NAMES IN PARAMSINORDER


% A script that will simulate fMRI BOLD data and fit a model with or
% without Q+ control
%
% Syntax:
%  [psiParamsFit]=simulate(model, paramsDomain, qpPres, varargin)
%
% Description:
%	Takes in a model and a possible set of parameters and whether or not Q+ 
%   is in control of things flag. 
%
% Inputs:
%   model                 - A function handle.
%                             @doeTemporalModel
%                             @watsonTemporalModel
%   paramsDomain          - Struct consisting of upper bounds, lower
%                           bounds, and intervals for all necessary parameters 
%                             DoE    (n=5): Sr, k1, k2, beta, sigma
%                             Watson (n=5): tau, kappa, zeta, beta, sigma
%
%   qpPres                - Logical: 
%                           true  - run simulation with Q+ stimulus presentation
%                           false - run simulation without Q+ stimulus presentation
%
% Optional key/value pairs (used in fitting):
%  
%p.addParameter('simulatedPsiParams',@isstruct);
%p.addParameter('headroom', 0.1, @isnumeric);
%p.addParameter('maxBOLD', 1.0, @isscalar);
%p.addParameter('maxBOLDSimulated', 1.5, @isscalar);
%p.addParameter('rngSeed',rng(1,'twister'),@isnumeric);
%p.addParameter('TR',800, @isnumeric);
%p.addParameter('trialLength',12, @isnumeric);
%p.addParameter('outNum','test',@ischar);
%p.addParameter('seed','choose');
%p.addParameter('nTrials',10,@isnumeric);
%p.addParameter('stimulusStructDeltaT',100,@isnumeric);
%p.addParameter('baselineStimulus',0);
%p.addParameter('maxBOLDStimulus',15);
%p.addParameter('nOutcomes',51,@isnumeric);
%p.addParameter('headroom',.1,@isscalar);
%
% Outputs:
%   psiParamsFit          - 1xn vector returning the BADS best fit for the
%                           parameters
%   maxBOLD               - Scalar. Best estimate at the maximum BOLD
%                           value.


%Example: 
%{
model = @doeTemporalModel;

paramsDomain = struct;
paramsDomain.Sr = 0.899:0.025:1.099;
paramsDomain.k1 = 0.01:0.04:0.4;
paramsDomain.k2 = 0.01:0.04:0.4;
paramsDomain.beta = 0.8:0.1:1.4; % Amplitude of the scaled response; should converge to unity
paramsDomain.sigma = 0.3:0.2:1;	% Standard deviation of the scaled (0-1) noise

simulatedPsiParams = struct;
simulatedPsiParams.Sr = .98;
simulatedPsiParams.k1 = .04;
simulatedPsiParams.k2 = .06;
simulatedPsiParams.beta = 1;
simulatedPsiParams.sigma = .1;

qpPres = logical(0);

[psiParamsFit]=simulate(model, paramsDomain, qpPres,...,
 'simulatedPsiParams', simulatedPsiParams);

%}

p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));
p.addRequired('paramsDomain',@isstruct);
p.addRequired('qpPres',@islogical);

% Optional params
p.addParameter('simulatedPsiParams',@isstruct);
p.addParameter('headroom', 0.1, @isnumeric);
p.addParameter('maxBOLD', 1.0, @isscalar);
p.addParameter('maxBOLDSimulated', 1.5, @isscalar);
p.addParameter('rngSeed',rng(1,'twister'),@isnumeric);
p.addParameter('TR',800, @isnumeric);
p.addParameter('trialLength',12, @isnumeric);
p.addParameter('outNum','test',@ischar);
p.addParameter('seed','choose');
p.addParameter('nTrials',10,@isnumeric);
p.addParameter('stimulusStructDeltaT',100,@isnumeric);
p.addParameter('baselineStimulus',0);
p.addParameter('maxBOLDStimulus',15);
p.addParameter('nOutcomes',51,@isnumeric);
p.addParameter('noiseSD',.1,@isscalar);
p.addParameter('stimulusDomain',@iscell);

% Parse
p.parse( model, paramsDomain, qpPres, varargin{:});

% Some variable name cleaning
headroom = p.Results.headroom;
maxBOLD = p.Results.maxBOLD;
maxBOLDSimulated = p.Results.maxBOLDSimulated;
trialLength = p.Results.trialLength;
TR = p.Results.TR;
outNum = p.Results.outNum;
seed = p.Results.seed;
nTrials = p.Results.nTrials; 
stimulusStructDeltaT = p.Results.stimulusStructDeltaT; 
baselineStimulus = p.Results.baselineStimulus;
maxBOLDStimulus = p.Results.maxBOLDStimulus;
noiseSD = p.Results.noiseSD;



% Establish qpParams
myQpParams = qpParams;


% Model specific processing. This should be functionalized. 
modelsCreated = {'doeTemporalModel','watsonTemporalModel'};
s = functions(model);
assert(any(strcmp(modelsCreated,s.function)==1),'Model not defined');

if contains(s.function,'doe')
    paramNamesInOrder = {'Sr', 'k1', 'k2', 'beta', 'sigma'};
    for i = 1:length(paramNamesInOrder)
        assert(isfield(paramsDomain,paramNamesInOrder{i}),'Parameter missing or misnamed for DoE model.');
    end
    % Create an anonymous function from qpDoETemporalModel in which we
    % specify the number of outcomes for the y-axis response
    myQpParams.qpPF = @(f,p) qpDoETemporalModel(f,p,myQpParams.nOutcomes,headroom);
    
elseif contains(s.function,'watson')
    paramNamesInOrder = {'tau', 'kappa', 'zeta', 'beta', 'sigma'};
    for i = 1:length(paramNamesInOrder)
        assert(isfield(paramsDomain,paramNamesInOrder{i}),'Parameter missing or misnamed for Watson model.');
    end
    % Create an anonymous function from qpDoETemporalModel in which we
    % specify the number of outcomes for the y-axis response
    myQpParams.qpPF = @(f,p) qpWatsonTemporalModel(f,p,myQpParams.nOutcomes,headroom);
    
end



% If seed is a keyword, pick a random seed.
if strcmp(seed,'choose')
    rng('shuffle'); seed = randi(2^32);
end


% Pick some random params to simulate if none provided (but set the neural
% noise to .1 SD and beta = 1)
if isempty(p.Results.simulatedPsiParams)
    simulatedPsiParams = struct;
    for i = 1:length(paramNamesInOrder)
        simulatedPsiParams.(paramNamesInOrder{i}) = randsample(paramsDomain.(paramNamesInOrder{i}),1);
    end
    % Beta is always one
    simulatedPsiParams.beta = 1;
    simulatedPsiParams.sigma = noiseSD;
else
    simulatedPsiParams = p.Results.simulatedPsiParams;
    % Beta will converge to 1 as maxBOLD gets closer and closer to the
    % simulated maxBOLD. As a result, when simulating data, beta should always
    % be set to 1. And, Q+ should always be able to incorporate 1 in its
    % domain. Assert these conditions are true. 
    assert(simulatedPsiParams.beta==1,'Simulated Beta should always be 1.');
    assert(ismember(1,paramsDomain.beta),'The domain for beta should always include 1.');
end


%% Some initialization
% Add the stimulus domain. 
if isempty(p.Results.stimulusDomain)
    %~Log spaced frequencies between 0 and 30 Hz
    myQpParams.stimParamsDomainList = {[baselineStimulus,1.875,3.75,7.5,15,30,60]};
end

% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = zeros(1,length(paramNamesInOrder));
upperBounds = zeros(1,length(paramNamesInOrder));
for i = 1:length(paramNamesInOrder)
    lowerBounds(i) = paramsDomain.(paramNamesInOrder{i})(1);
    upperBounds(i) = paramsDomain.(paramNamesInOrder{i})(end);
end

% Constrain bounds on beta to be very tight around 1.
lowerBounds{find(strcmp(paramNamesInOrder,'beta' ))} = .999;
upperBounds{find(strcmp(paramNamesInOrder,'beta' ))} = 1.001;



% We also want to make sure that the veridical values are actually within
% the domain bounds.
for param = 1:length(simulatedPsiParams)
    assert(lowerBounds(param) < simulatedPsiParams(param) && upperBounds(param) > simulatedPsiParams(param),...
        'Parameter %d is not within the bounds of the parameter domain.',param);
end

% Create and save an rng seed to use for this simulation.
rngSeed = rng(seed);
rngSeed = rng(seed);


% Create a simulated observer with binned output
myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,simulatedPsiParams);

% Initialize Q+
questData = qpInitialize(myQpParams);

% Tack on a continuous output simulated observer to myQpParams
myQpParams.continuousPF = @(f) model(f,simulatedPsiParams);

% Create a copy of Q+
questDataUntrained = questData;

% Create a stimulusVec to hold the trial across the loops
stimulusVec = nan(1,nTrials);

%% Run simulated trials
for tt = 1:nTrials
    
    % If it is the first three trials we force a baseline or maxBOLD event
    if tt<=3
        if tt == 1 || tt == 3
            stimulusVec(tt) = baselineStimulus;
            fprintf('Initial baseline stimulus: %0.3f\n',stimulusVec(tt));
        else
            stimulusVec(tt) = maxBOLDStimulus;
            fprintf('Initial maxBOLD stimulus: %0.3f\n',stimulusVec(tt));
        end
    else
        if ~simulateConstantStimuli
            % get random stimulus
            stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            fprintf('Stimuli chosen randomly: %0.3f\n',stimulusVec(tt));
        else
            % get next stimulus from Q+
            stimulusVec(tt) = qpQuery(questData);
            fprintf('Stimuli chosen by Q+: %0.3f\n',stimulusVec(tt));
        end
    end
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % that could be evoked by a stimulus (relative to the baseline
    % stimulus), which is the beta value of the model
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if tt > 2
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        maxBOLD = maxBOLD.*psiParamsQuest(4);
        fprintf('Using Q+ fit to generate maxBOLD \nmaxBOLD = %0.3f\n Q+ parameters: %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
        maxBOLD, psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));
    end

    % Create a packet
    thePacket = createPacket('nTrials',tt,...,
        'trialLengthSecs',trialLength,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);

    % Obtain outcomes from tfeUpdate 
    [outcomes] = ...
        tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
        'maxBOLDSimulated',maxBOLDSimulated,...
        'rngSeed',rngSeed.Seed,...,
        'maxBOLD',maxBOLD,...,
        'TRmsecs', TR,...,
        'noiseSD',simulatedPsiParams(5));

    % Grab a naive copy of questData
    questData = questDataUntrained;

    % Update quest data structure. This is the slow step in the simulation.
    for yy = 1:tt
        questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
    end
       
    
end

%% Adjust Beta/MaxBOLD tradeoff
% For the final parameter estimate, we want to assume that Beta is 1 and
% maxBOLD is whatever it WOULD BE if beta were 1. 

% Grab our current beta estimate is: 
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
betaGuess = psiParamsQuest(4);

% Divide maxBOLD by our beta estimate: (beta / beta) = 1, so
% new maxBOLD = maxBOLD/beta 
maxBOLD = maxBOLD.*betaGuess;

% Now run through the fitting steps again with the new maxBOLD
% Create a packet
thePacket = createPacket('nTrials',tt,...,
    'trialLengthSecs',trialLength,...,
    'stimulusStructDeltaT',stimulusStructDeltaT);

% Obtain outcomes from tfeUpdate
[outcomes] = ...
    tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
    'maxBOLDSimulated',maxBOLDSimulated,...
    'rngSeed',rngSeed.Seed,...,
    'maxBOLD',maxBOLD,...,
    'TRmsecs', TR,...,
    'noiseSD',simulatedPsiParams(5));

% Grab a naive copy of questData
questData = questDataUntrained;

% Update quest data structure. This is the slow step in the simulation.
for yy = 1:tt
    questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
end



%% Print some final output to the log

% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
fprintf('Simulated parameters:              %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
    simulatedPsiParams(1),simulatedPsiParams(2),simulatedPsiParams(3),simulatedPsiParams(4),simulatedPsiParams(5));
fprintf('FINAL Max posterior QUEST+ parameters:   %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));

psiParamsBads = psiParamsQuest;
psiParamsBads(4) = 1;

% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsBads,questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds,...
    'plausibleLowerBounds',lowerBounds,'plausibleUpperBounds',upperBounds);
fprintf('FINAL Maximum likelihood fit parameters: %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4),psiParamsFit(5));
fprintf('FINAL maxBOLD estimate: %0.3f',maxBOLD);

%% Output
T = array2table(psiParamsFit,'VariableNames',paramNamesInOrder);
T.maxBOLD = maxBOLD;
outfilename = horzcat(outNum,'.csv');
%save(outfilename,psiParamsFit);
writetable(T,outfilename);

end




