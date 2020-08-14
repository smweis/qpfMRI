function [psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain, varargin)
%% [psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain, varargin)
% A script that will simulate fMRI BOLD data and fit a model with or
% without Q+ control
%
% Syntax:
%  [psiParamsFit]=simulate(model, paramsDomain, varargin)
%
% Description:
%	Takes in a model and a possible set of parameters and whether or not Q+ 
%   is in control of things flag. 
%
% Inputs:
%   model                 - A function handle. This should be the
%                           'continuous function.' The Quest+ specific
%                           function will be defined in the model-specific
%                           code block. Currently supported:
%                             @doeTemporalModel
%                             @watsonTemporalModel
%   paramsDomain          - Struct consisting of upper bounds, lower
%                           bounds, and intervals for all necessary
%                           parameters. All models should have beta and
%                           sigma as parameters.
%                             DoE    (n=5): Sr, k1, k2, beta, sigma
%                             Watson (n=5): tau, kappa, zeta, beta, sigma
%
% Optional key/value pairs:
%
%   'qpPres'                -  Logical: (Default = false)
%                              true  - run simulation with Q+ stimulus choice
%                              false - run simulation with random stimulus
%                                      choice.
%   'simulatedPsiParams'    - Struct: (Default = randomly selected values
%                             from paramsDomain). The veridical parameters
%                             that are used for the forward model. Beta must
%                             be 1. 
%	'headroom'              - Scalar: (Default = 0.1)
%                             The proportion of the nOutcomes from qpParams 
%                             that will be used as extra on top and bottom.
%   'maxBOLD'               - Scalar: (Default = 1.0)
%                             The initial guess for the maximum BOLD value
%                             with respect to the baseline stimulus. 
%   'maxBOLDSimulated'      - Scalar: (Default = 1.5)
%                             The value (in % change units) of the
%                             maximum expected response to a stimulus w.r.t.
%                             the response to the baseline stimulus.
%	'seed'                  - No check (Default = 'choose')
%                             The value to initialize the rng. If 'choose'
%                             it will randomly initialize the seed using
%                             'shuffle'. 
%	'TR'                    - Integer: (Default = 800)
%                             Length of the time to repetition (TR) in
%                             milliseconds. 
%	'trialLength'           - Integer: (Default = 12)
%                             Length of one trial in seconds.
%	'outNum'                - String: (Default = 'test')
%                             Name of the output file (e.g., 'test.csv')
%	'outFolder'             - String: (Default = 'Results')
%                             Name of the output file (e.g., './Results')
%	'nTrials'               - Integer (Default = 10) 
%                             Number of trials to simulate. 
%	'stimulusStructDeltaT'  - Integer (Default = 100)
%                             Resolution of the stimulus struct in milliseconds
%                             (e.g., a value will be created every 100 ms).
%	'baselineStimulus'      - No check (Default = 0)
%                             tfeUpdate requires a baseline stimulus for
%                             which every value will be referenced. 
%	'maxBOLDStimulus'       - No check (Default = 15)
%                             It may help maxBOLD to require an early trial
%                             be a value we expect would result in a large
%                             BOLD response. 
%	'nOutcomes'             - Integer (Default = 51)
%                             The number of outcome bins Q+ can assign a
%                             response to. The larger this is the slower
%                             Q+ will be. 
%   'noiseSD'               - 1xn vector or scalar (Default = .1)
%                             Must be relative to maxBOLDSimulated. 
%                             In the absence of a simulatedPsiParam.sigma,
%                             noiseSD will allow the specification of a
%                             specific noise value (or selection from a
%                             vector of values). 
%	'stimulusDomain'        - Cell array (Default = a range of stimulus values).
%                             All possible stimulus values that can be
%                             assigned.
%   'stimulusDomainSpacing  - Default = 'lin'. 
%                             Whether the stimulusDomain is spaced linear
%                             or log.
%	'questDataCopy'         - Struct (Default is empty)
%                             If it is not necessary to reinitialize Q+
%                             (questData untrained is in memory), you can
%                             pass the initialized questDataCopy to
%                             simulate to save the time initializing. Note,
%                             any change to paramsDomain requires
%                             re-initialization. Be cautious changes other
%                             options without re-initializing. 
% Optional key/value pairs (used in plotting):
%   'showPlots'             - Logical: (Default = false)
%                             Whether to show plots.
%   'figWidth'              - Integer (Default = 900)
%                             Width of figure window size.
%   'figHeight'             - Integer (Default = 900)
%                             Height of figure window size.
%   'saveGif'               - Logical (Default = false)
%                             Whether to save the animated plot as a gif.
% Outputs:
%   psiParamsFit          - 1xn vector returning the BADS best fit for the
%                           parameters
%   maxBOLD               - Scalar. Best estimate at the maximum BOLD
%                           value.
%   questDataCopy         - Struct. Copy of initialized questData.

%Example: 
%{
---------------------------------------------------------------------------
Example 1: Logistic Model

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
nTrials = 100;

% Allow Q+ to control the stimuli or not (false).
qpPres = true;

nOutcomes = 15;% Set the number of outcome categories / bins.

showPlots = true; % Do you want to see plots?


% The range of BOLD signal to simulate (e.g., from baseline to maximum BOLD)
maxBOLDSimulated = 1.5;
% How noisy simulated BOLD data are in units of maxBOLDSimulated
noiseSD = .02; 
%How long the trials are (in seconds).
trialLength = 12;

% Note, this will save a copy of questData after it is initialized. 
[psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain,...,
'qpPres',qpPres, 'showPlots',showPlots,'stimulusDomain',stimulusDomain,...,
'stimulusDomainSpacing',stimulusDomainSpacing,'noiseSD',noiseSD,'nTrials',nTrials,...,
'maxBOLDSimulated',maxBOLDSimulated,'trialLength',trialLength,...,
'nOutcomes',nOutcomes);
---------------------------------------------------------------------------
Example 1 Time Saver code:

% Time saver for debugging: After running one of the above examples, keep
% everything in memory and run the line below. Especially useful if the
% paramsDomain is large or multi-dimensional.

% If paramsDomain is not changed, the following line can be run with
% questDataCopy as an optional argument to save the initialization step.
[psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain,...,
'qpPres',qpPres, 'showPlots',showPlots,'stimulusDomain',stimulusDomain,...,
'stimulusDomainSpacing',stimulusDomainSpacing,'noiseSD',noiseSD,...,
'simulatedPsiParams',simulatedPsiParams,'nTrials',nTrials,...,
'maxBOLDSimulated',maxBOLDSimulated,'trialLength',trialLength,...,
'nOutcomes',nOutcomes,'questDataCopy',questDataCopy);



%}

%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));
p.addRequired('paramsDomain',@isstruct);

% Optional params
p.addParameter('qpPres',false,@islogical);
p.addParameter('simulatedPsiParams',{});
p.addParameter('headroom', 0.1, @isnumeric);
p.addParameter('maxBOLD', 1.0, @isscalar);
p.addParameter('maxBOLDSimulated', 1.5, @isscalar);
p.addParameter('TR',800, @isnumeric);
p.addParameter('trialLength',12, @isnumeric);
p.addParameter('outNum','test',@ischar);
p.addParameter('outFolder','Results',@ischar);
p.addParameter('seed','choose');
p.addParameter('nTrials',10,@isnumeric);
p.addParameter('stimulusStructDeltaT',100,@isnumeric);
p.addParameter('baselineStimulus','');
p.addParameter('maxBOLDStimulus','');
p.addParameter('nOutcomes',5,@isnumeric);
p.addParameter('noiseSD',.1,@isvector);
p.addParameter('stimulusDomain',{},@iscell);
p.addParameter('stimulusDomainSpacing','lin',@ischar);
p.addParameter('questDataCopy',{},@isstruct);

% Optional params for plotting
p.addParameter('showPlots',false,@islogical);
p.addParameter('figWidth',1100,@isnumeric);
p.addParameter('figHeight',1100,@isnumeric);
p.addParameter('saveGif',false,@islogical);

% Parse
p.parse( model, paramsDomain, varargin{:});


% Establish qpParams
myQpParams = qpParams();

% Some variable name cleaning
qpPres = p.Results.qpPres;
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
% Put noiseSD on the scale of maxBOLDSimulated
noiseSD = p.Results.noiseSD .* maxBOLDSimulated;
myQpParams.nOutcomes = p.Results.nOutcomes;

showPlots = p.Results.showPlots;
figWidth = p.Results.figWidth;
figHeight = p.Results.figHeight;

try
    questDataCopy = p.Results.questDataCopy;
catch
    warning('No questData found. Will initialize, which will take some time.');
end

%% Check the model is supported and correct and return the model-specific values. 
[paramNamesInOrder, myQpParams.qpPF, myQpParams.psiParamsDomainList] = checkModel(model,...,
    paramsDomain,'nOutcomes',myQpParams.nOutcomes,'headroom',headroom);

%% Handle inputs
% This finds beta and sigma no matter where they are.
betaIndex = find(strcmp(paramNamesInOrder,'beta'));
sigmaIndex = find(strcmp(paramNamesInOrder,'sigma'));

% If seed is a keyword, pick a random seed. We use shuffle so it's new
% every time.
if strcmp(seed,'choose')
    rngSeed = rng('shuffle'); rngSeed = rng('shuffle');
else 
    fprintf('Initial random call %.04f\n',rand);
    rngSeed = rng(seed); rngSeed = rng(seed);
    fprintf('Random call after seed %.04f\n',rand);
end

%% Add the stimulus domain.  
if iscell(p.Results.stimulusDomain)
    myQpParams.stimParamsDomainList = p.Results.stimulusDomain;
elseif isvector(p.Results.stimulusDomain)
    myQpParams.stimParamsDomainList = {p.Results.stimulusDomain};
else
    warning('No stimulus domain specified. Defaulting to values between 0.01 and 1.');
    myQpParams.stimParamsDomainList = {makeDomain(.01,1,25)};
end

% Create baseline stimulus and maxBOLD stimulus if not passed.
if isempty(p.Results.baselineStimulus)
    baselineStimulus = min(myQpParams.stimParamsDomainList{1});
else
    baselineStimulus = p.Results.baselineStimulus;
end

if isempty(p.Results.maxBOLDStimulus)
    maxBOLDStimulus = max(myQpParams.stimParamsDomainList{1});
else
    maxBOLDStimulus = p.Results.maxBOLDStimulus;
end

%% Select the veridical psychometric paramteers for the function.
% Pick some random params to simulate if none provided but set beta to 1.
% We require the simulated parameters to result in a baseline trial = 0.
if isempty(p.Results.simulatedPsiParams)
    simulatedPsiParams = zeros(1,length(paramNamesInOrder));
    stillSearching = true;
    while stillSearching
        for i = 1:length(paramNamesInOrder)
            simulatedPsiParams(i) = randsample(paramsDomain.(paramNamesInOrder{i}),1);
        end
        
        % Beta is always one for simulations
        simulatedPsiParams(betaIndex) = 1;
        
        % Simulated noise is selected from a random sample of noiseSD
        simulatedPsiParams(sigmaIndex) = randsample(noiseSD,1);
        
        if abs(model(baselineStimulus,simulatedPsiParams)) < simulatedPsiParams(betaIndex)/10000 && ...
                abs(model(maxBOLDStimulus,simulatedPsiParams)) < 1 && ...
                abs(model(maxBOLDStimulus,simulatedPsiParams)) > .99
            stillSearching = false;
        end
    end
    

else
    simulatedPsiParams = zeros(1,length(paramNamesInOrder));
    
    % A relic of naming simulatedPsiParams noise parameter as sigma.
    for i = 1:length(paramNamesInOrder)
        if i == sigmaIndex && ~isfield(p.Results.simulatedPsiParams,'sigma')
            simulatedPsiParams(sigmaIndex) = noiseSD;
        else
            simulatedPsiParams(i) =  p.Results.simulatedPsiParams.(paramNamesInOrder{i});
        end
    end
    
    % Beta will converge to 1 as maxBOLD gets closer and closer to the
    % simulated maxBOLD. As a result, when simulating data, beta should always
    % be set to 1. 
    simulatedPsiParams(betaIndex) = 1;
    assert(simulatedPsiParams(betaIndex)==1,'Simulated Beta should always be 1.');
    
    % Select simulated psychometric parameters whose range of outputs
    % are between 0 and 1 for the model being used. 
    if abs(model(baselineStimulus,simulatedPsiParams)) < simulatedPsiParams(betaIndex)/10000
        warning('Simulated psychometric parameters will result in minimum values below 0.\nMin possible value = %.02f',abs(model(baselineStimulus,simulatedPsiParams)));
    elseif abs(model(baselineStimulus,simulatedPsiParams)) > .01
        warning('Simulated psychometric parameters will result in minimum values greater than 0.\nMin possible value = %.02f',abs(model(baselineStimulus,simulatedPsiParams)));
    elseif abs(model(maxBOLDStimulus,simulatedPsiParams)) < .99
        warning('Simulated psychometric parameters will result in maximum BOLD values below 1.\nMax possible value = %.02f',model(maxBOLDStimulus,simulatedPsiParams));
    end
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
lowerBoundsConstrained = lowerBounds;
upperBoundsConstrained = upperBounds;
lowerBoundsConstrained(betaIndex) = .999;
upperBoundsConstrained(betaIndex) = 1.001;

% Make sure 1 is a member of the beta parameter domains. If it's not, add
% it in. 
if ~ismember(1,paramsDomain.beta)
    warning('The domain for beta should always include 1. Adding 1 to beta parameter domain. This may create non-linear spacing.');
    paramsDomain.beta = sort(paramsDomain.beta);
    paramsDomain.beta = [paramsDomain.beta(paramsDomain.beta<1) 1 paramsDomain.beta(paramsDomain.beta>1)];
end

% Veridical values should be within domain bounds, except for sigma.
for param = 1:length(paramNamesInOrder)
    if ~strcmp(paramNamesInOrder{param},'sigma')
        assert(lowerBounds(param) <= simulatedPsiParams(param)...
            && upperBounds(param) >= simulatedPsiParams(param),...,
            'Parameter %s is not within the bounds of the parameter domain.',paramNamesInOrder{param});
    end
end

% Create a simulated observer with binned output
myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,simulatedPsiParams);

%% Initialize Q+. 
% Save some time if Q+ has already been initialized with the same
% parameters and the user has passed it as an input. 
if isstruct(questDataCopy)
    questData = questDataCopy;
else
    % Warn the user that initializing has started and may take a minute.
    tic
    fprintf('Initializing Q+. This may take a minute...\n');
    questData = qpInitialize(myQpParams);
    questDataCopy = questData;
    toc
end

% Calculate the lower headroom bin offset. We'll use this for plotting.
nLower = max([1 round(headroom*myQpParams.nOutcomes)]);
nUpper = max([1 round(headroom*myQpParams.nOutcomes)]);
nMid = myQpParams.nOutcomes - nLower - nUpper;

% Tack on a continuous output simulated observer to myQpParams
myQpParams.continuousPF = @(f) model(f,simulatedPsiParams);

% Create a copy of Q+
questDataUntrained = questData;

% Create a stimulusVec to hold the trial across the loops
stimulusVec = nan(1,nTrials);

%% Plotting features
% Select the plotting type based on stimulus domain spacing.
if strcmpi(p.Results.stimulusDomainSpacing,'log')
    plotFunc = @semilogx;
else
    plotFunc = @plot;
end

% Create a plot in which we can track the model progress
if showPlots
    % Create a fine version of the stimulus space.
    if strcmpi(p.Results.stimulusDomainSpacing,'log')
        stimulusDomainFine = logspace(log(min(myQpParams.stimParamsDomainList{1})),...,
            log(max(myQpParams.stimParamsDomainList{1})),100);
    else
        stimulusDomainFine = linspace(min(myQpParams.stimParamsDomainList{1}),...,
            max(myQpParams.stimParamsDomainList{1}),100);
    end
    
    %% Plot the parameter domains
    [paramsFig] = plotParamsDomain(model, paramsDomain, stimulusDomainFine,...,
        'nOutcomes',myQpParams.nOutcomes,'headroom',headroom,...,
        'stimulusDomainSpacing',p.Results.stimulusDomainSpacing,...,
        'figHeight',figHeight,'figWidth',figWidth);
    set(gcf,'color','w');
    hold off;
    
    % Create an empty packet for plotting
    thePacket = createPacket('nTrials',nTrials,...,
        'trialLengthSecs',trialLength,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);
    
    % Initialize the main figure
    mainFig = figure('Position',[10 10 figWidth figHeight]);
    set(gcf,'color','w');
    hold on;
    
    %% Subplot (bottom panel): Simulated fMRI data and model response.
    subplot(3,4,[9 10 11 12])
    currentBOLDHandleData = plot(thePacket.stimulus.timebase./1000,zeros(size(thePacket.stimulus.timebase)),'-k');
    hold on
    currentBOLDHandleFit = plot(thePacket.stimulus.timebase./1000,zeros(size(thePacket.stimulus.timebase)),'-r');
    xlim([min(thePacket.stimulus.timebase./1000) max(thePacket.stimulus.timebase)./1000]);
    % Use yMax to set the range for the plot.
    yMax = maxBOLDSimulated + maxBOLDSimulated*simulatedPsiParams(sigmaIndex)*2;
    ylim([-yMax yMax]);
    xlabel('time [seconds]','FontSize',20);
    ylabel('BOLD fMRI [% change]','FontSize',20);
    title('Simulated BOLD fMRI data','FontSize',20);
    set(gca,'box','off');
    ax = gca;
    ax.FontSize = 15;
    drawnow
    
    %% Subplot (left top panel): Veridical model, individual trials, current model fit.
    subplot(3,4,[1 2 5 6])
    set(gca,'Box','off');
    % TO DO: MIGHT WANT TO FIX THIS DOWN THE ROAD??
    % Predicted relative response is adjusted for the baseline and scaled
    % to maxBOLD
    predictedRelativeResponse = (model(stimulusDomainFine,simulatedPsiParams) - ...
        model(baselineStimulus,simulatedPsiParams))*maxBOLDSimulated;  
    plotFunc(stimulusDomainFine,predictedRelativeResponse,'-k');
    ylim([-0.5 maxBOLDSimulated+1]);
    xlabel('Stimulus Values');
    ylabel('Relative response amplitude');
    title('Estimate of Model');
    hold on
    % Scatter plot of the stimulus values
    currentOutcomesHandle = scatter(nan,nan);
    lastOutcomeHandle = scatter(nan,nan);
    % Veridical model plot.
    currentTTFHandle = plot(stimulusDomainFine,model(stimulusDomainFine,simulatedPsiParams),'-k');
    linePlot = plot(nan);
    
    
    %% Subplot (right top panel): Annotations.
    subplot(3,4,[3 4]);
    title('Trial Information','FontSize',25);
    set(gca,'visible','off');
    drawnow
    
    %% Subplot (right middle panel): Entropy by trial.
    subplot(3,4,[7 8]);
    entropyAfterTrial = nan(1,nTrials);
    xlabel('Trial number');
    ylabel('Entropy');
    currentEntropyHandle = plot(1:nTrials,entropyAfterTrial,'*k');
    xlim([1 nTrials]);
    set(gca,'box','off');
    drawnow
end

%% Print the simulated psychometric parameters:
fprintf('Simulated parameters:\n'); 
for i = 1:length(simulatedPsiParams)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},simulatedPsiParams(i));
end
fprintf('\n'); 

%% Run simulated trials
for tt = 1:nTrials

    % If it is the first three trials we force a baseline or maxBOLD event
    % Require 3 each. 
    if tt<=6
        if mod(tt,2) > 0 % Baseline trial
            stimulusVec(tt) = baselineStimulus;
            trialString = sprintf('\nTrial %d: Initial baseline\n',tt);
            stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(tt));
        else % maxBOLD trial
            stimulusVec(tt) = maxBOLDStimulus;
            trialString = sprintf('\nTrial %d: Initial max BOLD\n',tt);
            stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(tt));
        end
    % After that, every 5th trial will be a baseline or a maxBOLD trial.    
    elseif mod(tt,5) == 0 % Baseline trial
        stimulusVec(tt) = baselineStimulus;
        trialString = sprintf('\nTrial %d: Initial baseline\n',tt);
        stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(tt));
    elseif mod(tt,10) == 0 % maxBOLD trial
        stimulusVec(tt) = maxBOLDStimulus;
        trialString = sprintf('\nTrial %d: Initial max BOLD\n',tt);
        stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(tt));
    % Every other trial will be selected randomly, or by Q+
    else
        if ~qpPres
            % get random stimulus
            trialString = sprintf('\nTrial %d: Random Stimulus Selection\n',tt);
            stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(tt));
        else
            % get next stimulus from Q+
            trialString = sprintf('\nTrial %d: Q+ Stimulus Selection\n',tt);
            stimulusVec(tt) = qpQuery(questData);
            stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(tt));
        end
    end
    fprintf(trialString);
    fprintf(stimString);
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % that could be evoked by a stimulus (relative to the baseline
    % stimulus), which is the beta value of the model
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if tt > 2
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        % Update maxBOLD by the current best guess for Beta.
        maxBOLD = maxBOLD.*psiParamsQuest(betaIndex);
    end

    % Create a packet
    thePacket = createPacket('nTrials',tt,...,
        'trialLengthSecs',trialLength,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);

    % Obtain outcomes from tfeUpdate 
    [outcomes, modelResponseStruct, thePacketOut, ~, baselineEstimate] = ...
        tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
        'maxBOLDSimulated',maxBOLDSimulated,...
        'rngSeed',rngSeed,...,
        'maxBOLD',maxBOLD,...,
        'TRmsecs', TR,...,
        'noiseSD',simulatedPsiParams(sigmaIndex),...,
        'headroom',headroom);

    % Grab a naive copy of questData
    questData = questDataUntrained;

    % Update quest data structure. This is the slow step in the simulation.
    for yy = 1:tt
        questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
    end
    
    % TO DO: MIGHT NEED TO FIX THIS AS WELL. 
    % Get the yValues out of the outcome bins.
    yVals = (outcomes - nLower - 1)./nMid;
    psiParamsIndex = qpListMaxArg(questData.posterior);
    psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
    yVals = yVals*psiParamsQuest(betaIndex)*maxBOLD;
    yValsStim = yVals + baselineEstimate;
    
    % Print results to the console.
    resultString = sprintf('Output value %.03f\n',yVals(end));
    fprintf(resultString);
    entropyString = sprintf('Entropy after last trial: %.03f\n',questData.entropyAfterTrial(end));
    fprintf(entropyString);
    fprintf('\nQ+ parameters\n');
    for i = 1:length(psiParamsQuest)
        fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsQuest(i));
    end
    fprintf('\nmaxBOLD estimate = %0.3f\n',maxBOLD);
    fprintf('\n');
    
    %% Plot the ongoing results
    % Update the plots
    if showPlots
        % Delete all previous annotations.
        delete(findall(gcf,'type','annotation'));
        %% Subplot (bottom panel): Simulated fMRI data and model response.
        subplot(3,4,[9 10 11 12]);
        delete(currentBOLDHandleData)
        delete(currentBOLDHandleFit)
        currentBOLDHandleData = plot(thePacketOut.response.timebase./1000,thePacketOut.response.values,'.k');
        currentBOLDHandleFit = plot(modelResponseStruct.timebase./1000,modelResponseStruct.values,'-r');
        delete(linePlot);
        % This plots each stimulus as a blue line, with y-value
        % corresponding to simulated BOLD (as estimated by the nOutcome 
        % assignment and maxBOLD).
        for m = 1:tt
            linePlot(m) = plot([(m-1)*trialLength m*trialLength],[yValsStim(m) yValsStim(m)],'Color','b','LineWidth',4);
            linePlot(m).Color(4) = .2;
        end
        linePlot(m).Color(4) = 1;
        drawnow
        
        %% Subplot (left top panel): Veridical model, individual trials, current model fit.
        subplot(3,4,[1 2 5 6])
        % Current guess at the TTF, along with stims and outcomes
        set(gca,'Box','off');
        stimulusVecPlot = stimulusVec;
        stimulusVecPlot(stimulusVecPlot==0)=min(myQpParams.stimParamsDomainList{1});
        delete(currentOutcomesHandle);
        delete(lastOutcomeHandle);
        currentOutcomesHandle = scatter(stimulusVecPlot(1:tt-1),yVals(1:end-1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
        lastOutcomeHandle = scatter(stimulusVecPlot(tt),yVals(end),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',1,'HandleVisibility','off');
        predictedQuestRelativeResponse = (model(stimulusDomainFine,psiParamsQuest) - ...
            model(baselineStimulus,psiParamsQuest))*psiParamsQuest(betaIndex)*maxBOLD;
        delete(currentTTFHandle)
        currentTTFHandle = plotFunc(stimulusDomainFine,predictedQuestRelativeResponse,'-r');
        legend('Veridical Model','Individual Trials','Best-fit Model','Location','northwest');
        drawnow
        
        
        %% Subplot (right top panel): Annotations.
        subplot(3,4,[3 4]);
        title('Trial Information','FontSize',25);
        set(gca,'visible','off');
        ax = gca;
        xLoc = ax.Position(1);
        yLoc = ax.Position(2);
        annotation('textbox',[xLoc yLoc .4 .2],'String',sprintf('%s%s%s%s',trialString,stimString,resultString,entropyString),'EdgeColor','none','FontSize',15);
        drawnow
        
        %% Subplot (right middle panel): Entropy by trial.
        subplot(3,4,[7 8]);
        delete(currentEntropyHandle)
        entropyAfterTrial(1:tt)=questData.entropyAfterTrial;
        plot(1:nTrials,entropyAfterTrial,'*k');
        xlim([1 nTrials]);
        ylim([0 nanmax(entropyAfterTrial)]);
        set(gca,'box','off');
        title('Model entropy by trial number');
        xlabel('Trial number');
        ylabel('Entropy');
        drawnow
       
        hold off;
        
        % Save this plot as a GIF?
        if p.Results.saveGif
            if tt == 1
                % Requires gif package from fileexchange.
                try
                    gif('qpSimulate.gif','DelayTime',1,'frame',gcf,'nodither');
                catch
                    warning('You tried to save a gif, but perhaps you are missing the gif package.')
                end
            else
                gif;
            end
        end
    end
    
end

%% Adjust Beta/MaxBOLD tradeoff
% For the final parameter estimate, we want to assume that Beta is 1 and
% maxBOLD is whatever it WOULD BE if beta were 1. 

% Grab our current beta estimate is: 
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds,...
    'plausibleLowerBounds',lowerBounds,'plausibleUpperBounds',upperBounds);
maxBOLD = maxBOLD.*psiParamsFit(betaIndex);

% Now run through the fitting steps again with the new maxBOLD
thePacket = createPacket('nTrials',tt,...,
    'trialLengthSecs',trialLength,...,
    'stimulusStructDeltaT',stimulusStructDeltaT);

[outcomes] = ...
    tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
    'maxBOLDSimulated',maxBOLDSimulated,...
    'rngSeed',rngSeed,...,
    'maxBOLD',maxBOLD,...,
    'TRmsecs', TR,...,
    'noiseSD',simulatedPsiParams(sigmaIndex));

questData = questDataUntrained;
for yy = 1:tt
    questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
end

%% Print some final output to the log
fprintf('--------------------------------------------------------------\n');
fprintf('FINAL VALUES\n');
fprintf('--------------------------------------------------------------\n');
% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
% Print simulated parameters.
fprintf('Simulated parameters:              '); 
for i = 1:length(simulatedPsiParams)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},simulatedPsiParams(i));
end
% Print Q+ Fit parameters.
fprintf('\nMax posterior QUEST+ parameters:   '); 
for i = 1:length(psiParamsQuest)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsQuest(i));
end
% Prepare for BADS fit.
psiParamsBads = psiParamsQuest;
psiParamsBads(betaIndex) = 1;

% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsBads,questData.nOutcomes,...
    'lowerBounds', lowerBoundsConstrained,'upperBounds',upperBoundsConstrained,...
    'plausibleLowerBounds',lowerBoundsConstrained,'plausibleUpperBounds',upperBoundsConstrained);
% Print BADS best fit parameters.
fprintf('\nMax BADS parameters:               '); 
for i = 1:length(psiParamsFit)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsFit(i));
end
% Print maxBOLD
fprintf('\nmaxBOLD estimate: %0.3f\n',maxBOLD);

%% Final plot fits
if showPlots
    finalFig = figure('Position', [10 10 figWidth figHeight]);
    hold on;
    % TO DO: Might want to fix this for plotting.
    predictedBADSRelativeResponse = (model(stimulusDomainFine,psiParamsFit) - ...
        model(baselineStimulus,psiParamsFit)).*maxBOLD;
    plotFunc(stimulusDomainFine,predictedRelativeResponse,'-k','LineWidth',6);
    plotFunc(stimulusDomainFine,predictedBADSRelativeResponse,'-','Color','#FA4515','LineWidth',6);
    if strcmpi(p.Results.stimulusDomainSpacing,'log')
        set(gca,'XScale', 'log');
    end
    xlabel('Contrast');
    ylabel('Normalized Predicted Response');
    legend('Veridical Model','Final Model Fit','Location','northwest');
end


%% Output
% Two files will be created for each run. 
% [outnum]Results.csv contains the BADS results of the simulation. 
% [outnum]Params.csv contains the simulated parameters.

% Save figures
if showPlots
    delete(findall(gcf,'type','annotation'));
    set(mainFig,'Units','Inches');
    pos = get(mainFig,'Position');
    set(mainFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(mainFig,'simulate2.pdf','-dpdf','-r0');
    
    set(paramsFig,'Units','Inches');
    pos = get(paramsFig,'Position');
    set(paramsFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(paramsFig,'paramsFigs2.pdf','-dpdf','-r0');
    
    set(finalFig,'Units','Inches');
    pos = get(finalFig,'Position');
    set(finalFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(finalFig,'finalFig2.pdf','-dpdf','-r0');
end

resultsOut = array2table(psiParamsFit,'VariableNames',paramNamesInOrder);
resultsOut.maxBOLD = maxBOLD;
resultsFileName = [outNum 'Results.csv'];
resultsFolderName = ['.' filesep p.Results.outFolder filesep 'results'];

if ~exist(resultsFolderName,'dir')
    mkdir(resultsFolderName);
end

writetable(resultsOut,[resultsFolderName filesep resultsFileName]);

simulatedParamsOut = array2table(simulatedPsiParams,'VariableNames',paramNamesInOrder);
simulatedParamsOut.maxBOLDSimulated = maxBOLDSimulated;
simulatedParamsOut.nOutcomes = p.Results.nOutcomes;

paramsFileName = [outNum 'Params.csv'];
paramsFolderName = ['.' filesep p.Results.outFolder filesep 'params'];

if ~exist(paramsFolderName,'dir')
    mkdir(paramsFolderName);
end

writetable(simulatedParamsOut,[paramsFolderName filesep paramsFileName]);

end




