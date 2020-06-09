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

% Outputs:
%   psiParamsFit          - 1xn vector returning the BADS best fit for the
%                           parameters
%   maxBOLD               - Scalar. Best estimate at the maximum BOLD
%                           value.
%   questDataCopy         - Struct. Copy of initialized questData.

%Example: 
%{
---------------------------------------------------------------------------
% Example 1: DoE Temporal Model - a somewhat poorly behaved model.
model = @doeTemporalModel;

paramsDomain = struct;
paramsDomain.Sr = 0.899:0.025:1.099;
paramsDomain.k1 = 0.01:0.04:0.4;
paramsDomain.k2 = 0.01:0.04:0.4;
paramsDomain.beta = 0.8:0.1:1.4;
paramsDomain.sigma = 0.3:0.2:1.0;	

simulatedPsiParams = struct;
simulatedPsiParams.Sr = 1.004; 
simulatedPsiParams.k1 = .016; 
simulatedPsiParams.k2 = .118; 
simulatedPsiParams.beta = 1.0; 
simulatedPsiParams.sigma = .1;

qpPres = false;

showPlots = true;

[psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain,...,
'qpPres',qpPres,'simulatedPsiParams', simulatedPsiParams,...,
 'showPlots',showPlots);

---------------------------------------------------------------------------
Example 2: Logistic Model

model = @logistic;

paramsDomain = struct;
paramsDomain.slope = makeDomain(-1.2,-.2,10,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.beta = makeDomain(.75,1.25,11,'spacing','zeno');
paramsDomain.sigma = makeDomain(.05,2,10);

stimulusDomain = {makeDomain(.01,1,25)};
stimulusDomainSpacing = 'lin';
nTrials = 30;
qpPres = false;

showPlots = true;

simulatedPsiParams = struct;
simulatedPsiParams.slope = .12;
simulatedPsiParams.semiSat = .49;
simulatedPsiParams.beta = 1.0;
simulatedPsiParams.sigma = .4;

trialLength = 12;

% Note, this will save a copy of questData after it is initialized. 
[psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain,...,
'qpPres',qpPres, 'showPlots',showPlots,'stimulusDomain',stimulusDomain,...,
'stimulusDomainSpacing',stimulusDomainSpacing,...,
'simulatedPsiParams',simulatedPsiParams,'nTrials',nTrials,'trialLength',trialLength);
---------------------------------------------------------------------------
Time saver for debugging: After running one of the above examples, keep
everything in memory and run the line below. Especially useful if the
paramsDomain is large or multi-dimensional.

% If paramsDomain is not changed, the following line can be run with
questDataCopy as an optional argument to save the initialization step.
[psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain,...,
'qpPres',qpPres, 'showPlots',showPlots,'stimulusDomain',stimulusDomain,...,
'stimulusDomainSpacing',stimulusDomainSpacing,...,
'questDataCopy',questDataCopy,'simulatedPsiParams',simulatedPsiParams,'nTrials',nTrials);



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
p.addParameter('nOutcomes',51,@isnumeric);
p.addParameter('noiseSD',.1,@isvector);
p.addParameter('stimulusDomain',{},@iscell);
p.addParameter('stimulusDomainSpacing','lin',@ischar);
p.addParameter('questDataCopy',{},@isstruct);

% Optional params for plotting
p.addParameter('showPlots',false,@islogical);
p.addParameter('figWidth',1100,@isnumeric);
p.addParameter('figHeight',1100,@isnumeric);

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

% Changing nOutcomes has unpredictable effects on the search for sigma
if p.Results.nOutcomes ~= 51
    warning('51 bins (nOutcomes) are recommended. If you change this, be sure to alter the parameter domain for sigma.');
end

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
    myQpParams.stimParamsDomainList = {[baselineStimulus,1.875,3.75,7.5,15,30,60]};
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
        
        % simulatedPsiParams is selected from a random sample of
        % p.Results.noiseSD
        simulatedPsiParams(sigmaIndex) = randsample(p.Results.noiseSD,1);
        
        if abs(model(baselineStimulus,simulatedPsiParams)) < simulatedPsiParams(betaIndex)/10000
            stillSearching = false;
        elseif abs(model(maxBOLDStimulus,simulatedPsiParams)) < 1
            stillSearching = false;
        end
    end
    

else
    simulatedPsiParams = zeros(1,length(paramNamesInOrder));
    for i = 1:length(paramNamesInOrder)
        simulatedPsiParams(i) =  p.Results.simulatedPsiParams.(paramNamesInOrder{i});
    end
    % Beta will converge to 1 as maxBOLD gets closer and closer to the
    % simulated maxBOLD. As a result, when simulating data, beta should always
    % be set to 1. And, Q+ should always be able to incorporate 1 in its
    % domain. Assert these conditions are true. 
    % Beta is always one for simulations
    simulatedPsiParams(betaIndex) = 1;
        
    assert(simulatedPsiParams(betaIndex)==1,'Simulated Beta should always be 1.');
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




% We also want to make sure that the veridical values are actually within
% the domain bounds. Don't do this for sigma. 
for param = 1:length(paramNamesInOrder)
    if ~strcmp(paramNamesInOrder{param},'sigma')
        assert(lowerBounds(param) <= simulatedPsiParams(param)...
            && upperBounds(param) >= simulatedPsiParams(param),...,
            'Parameter %s is not within the bounds of the parameter domain.',paramNamesInOrder{param});
    end
end


% Create a simulated observer with binned output
myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,simulatedPsiParams);



% Initialize Q+. Save some time if we're debugging
if isstruct(questDataCopy)
    questData = questDataCopy;
else
    % Warn the user that we are initializing
    tic
    fprintf('Initializing Q+. This may take a minute...\n');
    questData = qpInitialize(myQpParams);
    questDataCopy = questData;
    toc
end

    
% Calculate the lower headroom bin offset. We'll use this later
nLower = round(headroom*myQpParams.nOutcomes);
nUpper = round(headroom*myQpParams.nOutcomes);
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
    if strcmpi(p.Results.stimulusDomainSpacing,'log')
        stimulusDomainFine = logspace(log(min(myQpParams.stimParamsDomainList{1})),...,
            log(max(myQpParams.stimParamsDomainList{1})),100);
    else
        stimulusDomainFine = linspace(min(myQpParams.stimParamsDomainList{1}),...,
            max(myQpParams.stimParamsDomainList{1}),100);
    end
    % First, we'll plot the parameter domains
    [paramsFig] = plotParamsDomain(model, paramsDomain, stimulusDomainFine,...,
        'nOutcomes',myQpParams.nOutcomes,'headroom',headroom,...,
        'stimulusDomainSpacing',p.Results.stimulusDomainSpacing,...,
        'figHeight',figHeight,'figWidth',figWidth);
    set(gcf,'color','w');
    
    hold off;
    
    % Create the packet early for plotting
    thePacket = createPacket('nTrials',nTrials,...,
        'trialLengthSecs',trialLength,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);
    
    % Initialize the figure
    mainFig = figure('Position',[10 10 figWidth figHeight]);
    set(gcf,'color','w');
    hold on;
    
    
    
    % Set up the BOLD fMRI response and model fit
    subplot(3,4,[9 10 11 12])
    currentBOLDHandleData = plot(thePacket.stimulus.timebase./1000,zeros(size(thePacket.stimulus.timebase)),'-k');
    hold on
    currentBOLDHandleFit = plot(thePacket.stimulus.timebase./1000,zeros(size(thePacket.stimulus.timebase)),'-r');
    xlim([min(thePacket.stimulus.timebase./1000) max(thePacket.stimulus.timebase)./1000]);
    yMax = maxBOLDSimulated + maxBOLDSimulated*simulatedPsiParams(sigmaIndex)*2;
    ylim([-yMax yMax]);
    xlabel('time [seconds]','FontSize',20);
    ylabel('BOLD fMRI [% change]','FontSize',20);
    title('Simulated BOLD fMRI data','FontSize',20);
    set(gca,'box','off');
    ax = gca;
    ax.FontSize = 15;
    
    % Set up the TTF figure
    subplot(3,4,[1 2 5 6])
    % MIGHT WANT TO FIX THIS DOWN THE ROAD??
    predictedRelativeResponse = model(stimulusDomainFine,simulatedPsiParams) - ...
        model(baselineStimulus,simulatedPsiParams);
    % May need to scale the predictedRelativeResponse here to account for
    % the offset produced by subtraction of the baseline amplitude.    
    plotFunc(stimulusDomainFine,predictedRelativeResponse,'-k');
    ylim([-0.5 1.5]);
    xlabel('Stimulus Values');
    ylabel('Relative response amplitude');
    title('Estimate of Model');
    hold on
    currentOutcomesHandle = scatter(nan,nan);
    currentTTFHandle = plot(stimulusDomainFine,model(stimulusDomainFine,simulatedPsiParams),'-k');
  
    % Set up the entropy x trial figure
    h = subplot(3,4,[3 4 7 8]);
    entropyAfterTrial = nan(1,nTrials);
    currentEntropyHandle = plot(1:nTrials,entropyAfterTrial,'*k');
    xlim([1 nTrials]);
    set(gca,'box','off');
    title('Model entropy by trial number');
    xlabel('Trial number');
    ylabel('Entropy');
end




%% Print the simulated psychometric parameters:
fprintf('Simulated parameters:\n'); 
for i = 1:length(simulatedPsiParams)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},simulatedPsiParams(i));
end
fprintf('\n'); 

%% Run simulated trials
for tt = 1:nTrials
    fprintf('\nTrial %d\n',tt);
    % If it is the first three trials we force a baseline or maxBOLD event
    % Require 3 each. 
    if tt<=6
        if mod(tt,2) > 0
            stimulusVec(tt) = baselineStimulus;
            string = sprintf('Stimulus (Initial baseline): %0.3f\n',stimulusVec(tt));
        else
            stimulusVec(tt) = maxBOLDStimulus;
            string = sprintf('Stimulus (Initial maxBOLD): %0.3f\n',stimulusVec(tt));
        end
    % Require a maxBOLD and a baseline every X trials (SAVE FOR LATER)
    else
        if ~qpPres
            % get random stimulus
            stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            string = sprintf('Stimulus (Random): %0.3f\n',stimulusVec(tt));
        else
            % get next stimulus from Q+
            stimulusVec(tt) = qpQuery(questData);
            string = sprintf('Stimulus (Q+): %0.3f\n',stimulusVec(tt));
        end
    end
    fprintf(string);
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % that could be evoked by a stimulus (relative to the baseline
    % stimulus), which is the beta value of the model
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if tt > 2
        
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
       
        maxBOLD = maxBOLD.*psiParamsQuest(betaIndex);
       
        fprintf('\nQ+ parameters\n');
        for i = 1:length(psiParamsQuest)
            fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsQuest(i));
        end
        
        fprintf('\nmaxBOLD estimate = %0.3f\n',maxBOLD);
        fprintf('\n');
    end

    % Create a packet
    thePacket = createPacket('nTrials',tt,...,
        'trialLengthSecs',trialLength,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);

    % Obtain outcomes from tfeUpdate 
    [outcomes, modelResponseStruct, thePacketOut] = ...
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
    
    % Get the yValues out of the outcome bins.
    yVals = (outcomes - nLower - 1)./nMid;
    psiParamsIndex = qpListMaxArg(questData.posterior);
    psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
    yVals = (yVals*psiParamsQuest(betaIndex))*maxBOLD;
    
    %% Plot the ongoing results
    % Update the plots
    if showPlots
        % Delete all previous annotations.
        delete(findall(gcf,'type','annotation'));
        % Simulated BOLD fMRI time-series and fit
        subplot(3,4,[9 10 11 12]);
        delete(currentBOLDHandleData)
        delete(currentBOLDHandleFit)
        currentBOLDHandleData = plot(thePacketOut.response.timebase./1000,thePacketOut.response.values,'.k');
        currentBOLDHandleFit = plot(modelResponseStruct.timebase./1000,modelResponseStruct.values,'-r');
        
        linePlot = plot([(tt-1)*trialLength tt*trialLength],[yVals(tt) yVals(tt)],'Color','b','LineWidth',4);
        linePlot.Color(4) = .2;
        drawnow
        
        
        % TTF figure
        subplot(3,4,[1 2 5 6])
        % Current guess at the TTF, along with stims and outcomes

        stimulusVecPlot = stimulusVec;
        stimulusVecPlot(stimulusVecPlot==0)=min(myQpParams.stimParamsDomainList{1});
        delete(currentOutcomesHandle);
        currentOutcomesHandle = scatter(stimulusVecPlot(1:tt),yVals,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
       
        predictedQuestRelativeResponse = model(stimulusDomainFine,psiParamsQuest) - ...
            model(baselineStimulus,psiParamsQuest);
        ax = gca;
        xLoc = ax.Position(1);
        yLoc = ax.Position(2);
        annotation('textbox',[xLoc+.6 yLoc-.3 .4 .2],'String',sprintf('Trial #%d\n%s',tt,string),'EdgeColor','none','FontSize',15);
        % MAY WANT TO FIX TO ALLOW LOG PLOTS
        delete(currentTTFHandle)
        currentTTFHandle = plot(stimulusDomainFine,predictedQuestRelativeResponse,'-r');
        legend('Veridical model','Stimulus outcome bins','Best-fit model','Location','northwest');
        drawnow
        
        
        % Entropy by trial
        h = subplot(3,4,[3 4 7 8]);
        delete(currentEntropyHandle)
        entropyAfterTrial(1:tt)=questData.entropyAfterTrial;
        plot(1:nTrials,entropyAfterTrial,'*k');
        xlim([1 nTrials]);
        ylim([0 nanmax(entropyAfterTrial)]);
        xlabel('Trial number');
        ylabel('Entropy');
        drawnow
        
        
        % Requires gif package from fileexchange.
        if tt == 1
            % Requires gif package from fileexchange.
            gif('simulate.gif','DelayTime',1,'frame',gcf,'nodither');
        else
            gif;
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
fprintf('Simulated parameters:              '); 
for i = 1:length(simulatedPsiParams)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},simulatedPsiParams(i));
end

fprintf('\nMax posterior QUEST+ parameters:   '); 

for i = 1:length(psiParamsQuest)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsQuest(i));
end


psiParamsBads = psiParamsQuest;
psiParamsBads(betaIndex) = 1;

% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsBads,questData.nOutcomes,...
    'lowerBounds', lowerBoundsConstrained,'upperBounds',upperBoundsConstrained,...
    'plausibleLowerBounds',lowerBoundsConstrained,'plausibleUpperBounds',upperBoundsConstrained);


fprintf('\nMax BADS parameters:               '); 
for i = 1:length(psiParamsFit)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsFit(i));
end

fprintf('\nmaxBOLD estimate: %0.3f\n',maxBOLD);

%% Final plot fits
if showPlots
    finalFig = figure('Position', [10 10 figWidth figHeight]);
    hold on;
    predictedBADSRelativeResponse = model(stimulusDomainFine,psiParamsFit) - ...
        model(baselineStimulus,psiParamsFit);
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




