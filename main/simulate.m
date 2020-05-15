function [psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain, qpPres, varargin)
%% simulate
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
%   qpPres                - Logical: 
%                           true  - run simulation with Q+ stimulus choice
%                           false - run simulation with random stimulus
%                                   choice.
%
% Optional key/value pairs (used in fitting):
%  %% TODO - I HAVE NOT DEFINED THESE YET!
%   'simulatedPsiParams',@isstruct);
%	'headroom', 0.1, @isnumeric);
%   'maxBOLD', 1.0, @isscalar);
%   'maxBOLDSimulated', 1.5, @isscalar);
%	'rngSeed',rng(1,'twister'),@isnumeric);
%	'TR',800, @isnumeric);
%	'trialLength',12, @isnumeric);
%	'outNum','test',@ischar);
%	'seed','choose');
%	'nTrials',10,@isnumeric);
%	'stimulusStructDeltaT',100,@isnumeric);
%	'baselineStimulus',0);
%	'maxBOLDStimulus',15);
%	'nOutcomes',51,@isnumeric);
%	'headroom',.1,@isscalar);
%	'stimulusDomain',{},@iscell);
%	'questDataCopy',{},@isstruct);
% Optional key/value pairs (used in plotting):
%   'showPlots',false,@islogical);
%   'minStim',.01,@isscalar);
%   'maxStim',100,@isscalar);
%   'figWidth',900,@isnumeric);
%   'figHeight',900,@isnumeric);

% Outputs:
%   psiParamsFit          - 1xn vector returning the BADS best fit for the
%                           parameters
%   maxBOLD               - Scalar. Best estimate at the maximum BOLD
%                           value.
% 

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

qpPres = false;

showPlots = true;

% Note, this will save a copy of questData after it is initialized. 
[psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain, qpPres,...,
 'simulatedPsiParams', simulatedPsiParams, 'showPlots',showPlots);

% If paramsDomain is not changed, the following line can be run with
questDataCopy as an optional argument to save the initialization step.
[psiParamsFit,maxBOLD,questDataCopy]=simulate(model, paramsDomain, qpPres,...,
 'simulatedPsiParams', simulatedPsiParams,'questDataCopy',questDataCopy);


%}

%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));
p.addRequired('paramsDomain',@isstruct);
p.addRequired('qpPres',@islogical);

% Optional params
p.addParameter('simulatedPsiParams',{},@isstruct);
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
p.addParameter('stimulusDomain',{},@iscell);
p.addParameter('questDataCopy',{},@isstruct);

% Optional params for plotting
p.addParameter('showPlots',false,@islogical);
p.addParameter('minStim',.01,@isscalar);
p.addParameter('maxStim',100,@isscalar);
p.addParameter('figWidth',900,@isnumeric);
p.addParameter('figHeight',900,@isnumeric);

% Parse
p.parse( model, paramsDomain, qpPres, varargin{:});


% Establish qpParams
myQpParams = qpParams();

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
myQpParams.nOutcomes = p.Results.nOutcomes;
showPlots = p.Results.showPlots;
minStim = p.Results.minStim;
maxStim = p.Results.maxStim;
figWidth = p.Results.figWidth;
figHeight = p.Results.figHeight;

try
    questDataCopy = p.Results.questDataCopy;
catch
    warning('No questData found. Will initialize, which will take some time.');
end

%% This function contains all supported models and returns the model-specific values. 
[paramNamesInOrder, myQpParams.qpPF, myQpParams.psiParamsDomainList] = checkModel(model,paramsDomain,myQpParams,headroom);

%% Handle more inputs
% This finds beta and sigma no matter where it is.
betaIndex = find(strcmp(paramNamesInOrder,'beta'));
sigmaIndex = find(strcmp(paramNamesInOrder,'sigma'));

% If seed is a keyword, pick a random seed.
if strcmp(seed,'choose')
    rng('shuffle'); seed = randi(2^32);
end

% Pick some random params to simulate if none provided (but set the neural
% noise to .1 SD and beta = 1)
if isempty(p.Results.simulatedPsiParams)
    simulatedPsiParams = zeros(1,length(paramNamesInOrder));
    for i = 1:length(paramNamesInOrder)
        simulatedPsiParams(i) = randsample(paramsDomain.(paramNamesInOrder{i}),1);
    end
    % Beta is always one
    simulatedPsiParams(betaIndex) = 1;
    simulatedPsiParams(sigmaIndex) = noiseSD;
else
    simulatedPsiParams = zeros(1,length(paramNamesInOrder));
    for i = 1:length(paramNamesInOrder)
        simulatedPsiParams(i) =  p.Results.simulatedPsiParams.(paramNamesInOrder{i});
    end
    % Beta will converge to 1 as maxBOLD gets closer and closer to the
    % simulated maxBOLD. As a result, when simulating data, beta should always
    % be set to 1. And, Q+ should always be able to incorporate 1 in its
    % domain. Assert these conditions are true. 
    assert(simulatedPsiParams(betaIndex)==1,'Simulated Beta should always be 1.');
    assert(ismember(1,paramsDomain.beta),'The domain for beta should always include 1.');
end

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
lowerBounds(betaIndex) = .999;
upperBounds(betaIndex) = 1.001;


% We also want to make sure that the veridical values are actually within
% the domain bounds. Don't do this for sigma. 
for param = 1:length(paramNamesInOrder)
    if ~strcmp(paramNamesInOrder{param},'sigma')
        assert(lowerBounds(param) < simulatedPsiParams(param)...
            && upperBounds(param) > simulatedPsiParams(param),...,
            'Parameter %s is not within the bounds of the parameter domain.',paramNamesInOrder{param});
    end
end

% Create and save an rng seed to use for this simulation.
rngSeed = rng(seed);
rngSeed = rng(seed);


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



% Tack on a continuous output simulated observer to myQpParams
myQpParams.continuousPF = @(f) model(f,simulatedPsiParams);

% Create a copy of Q+
questDataUntrained = questData;

% Create a stimulusVec to hold the trial across the loops
stimulusVec = nan(1,nTrials);

%% Plotting features

% Create a plot in which we can track the model progress
if showPlots
    
    
    % First, we'll plot the parameter domains
    plotParamsDomain(model, paramsDomain, 'minStim',minStim,'maxStim',maxStim,...,
        'figHeight',figHeight,'figWidth',figWidth);
    hold off;
    
    % Create the packet early
    thePacket = createPacket('nTrials',nTrials,...,
        'trialLengthSecs',trialLength,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);

    figure('Position',[10 10 figWidth figHeight]);
    hold on;
    % Set up the BOLD fMRI response and model fit
    subplot(3,1,1)
    currentBOLDHandleData = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-k');
    hold on
    currentBOLDHandleFit = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-r');
    xlim([min(thePacket.stimulus.timebase) max(thePacket.stimulus.timebase)]);
    ylim([-2 2]);
    xlabel('time [msecs]');
    ylabel('BOLD fMRI % change');
    title('BOLD fMRI data');
    
    % Set up the TTF figure
    subplot(3,1,2)
    freqDomain = logspace(log10(minStim),log10(maxStim),100);
    predictedRelativeResponse = model(freqDomain,simulatedPsiParams) - ...
        model(baselineStimulus,simulatedPsiParams);
    % May need to scale the predictedRelativeResponse here to account for
    % the offset produced by subtraction of the baseline amplitude.    
    semilogx(freqDomain,predictedRelativeResponse,'-k');
    ylim([-0.5 1.5]);
    xlabel('log stimulus Frequency [Hz]');
    ylabel('Relative response amplitude');
    title('Estimate of Model TTF');
    hold on
    currentOutcomesHandle = scatter(nan,nan);
    currentTTFHandle = plot(freqDomain,model(freqDomain,simulatedPsiParams),'-k');
    
    % Calculate the lower headroom bin offset. We'll use this later
    nLower = round(headroom*myQpParams.nOutcomes);
    nUpper = round(headroom*myQpParams.nOutcomes);
    nMid = myQpParams.nOutcomes - nLower - nUpper;
        
    % Set up the entropy x trial figure
    subplot(3,1,3)
    entropyAfterTrial = nan(1,nTrials);
    currentEntropyHandle = plot(1:nTrials,entropyAfterTrial,'*k');
    xlim([1 nTrials]);
    title('Model entropy by trial number');
    xlabel('Trial number');
    ylabel('Entropy');
end


%% Run simulated trials
for tt = 1:nTrials
    fprintf('\nTrial %d\n',tt);
    % If it is the first three trials we force a baseline or maxBOLD event
    if tt<=3
        if tt == 1 || tt == 3
            stimulusVec(tt) = baselineStimulus;
            fprintf('STIMULUS (Initial baseline): %0.3f\n',stimulusVec(tt));
        else
            stimulusVec(tt) = maxBOLDStimulus;
            fprintf('STIMULUS (Initial maxBOLD): %0.3f\n',stimulusVec(tt));
        end
    else
        if ~qpPres
            % get random stimulus
            stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            fprintf('STIMULUS (Random): %0.3f\n',stimulusVec(tt));
        else
            % get next stimulus from Q+
            stimulusVec(tt) = qpQuery(questData);
            fprintf('STIMULUS (Q+): %0.3f\n',stimulusVec(tt));
        end
    end
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % that could be evoked by a stimulus (relative to the baseline
    % stimulus), which is the beta value of the model
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if tt > 2
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        maxBOLD = maxBOLD.*psiParamsQuest(betaIndex);
        fprintf('maxBOLD estimate = %0.3f\n',maxBOLD);
        for i = 1:length(psiParamsQuest)
            fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsQuest(i));
        end
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
        'rngSeed',rngSeed.Seed,...,
        'maxBOLD',maxBOLD,...,
        'TRmsecs', TR,...,
        'noiseSD',simulatedPsiParams(sigmaIndex));

    % Grab a naive copy of questData
    questData = questDataUntrained;

    % Update quest data structure. This is the slow step in the simulation.
    for yy = 1:tt
        questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
    end
        % Update the plots
    if showPlots
        
        % Simulated BOLD fMRI time-series and fit
        subplot(3,1,1)
        delete(currentBOLDHandleData)
        delete(currentBOLDHandleFit)
        currentBOLDHandleData = plot(thePacketOut.response.timebase,thePacketOut.response.values,'.k');
        currentBOLDHandleFit = plot(modelResponseStruct.timebase,modelResponseStruct.values,'-r');
                
        % TTF figure
        subplot(3,1,2)
        % Current guess at the TTF, along with stims and outcomes
        yVals = (outcomes - nLower - 1)./nMid;
        stimulusVecPlot = stimulusVec;
        stimulusVecPlot(stimulusVecPlot==0)=minStim;
        delete(currentOutcomesHandle);
        currentOutcomesHandle = scatter(stimulusVecPlot(1:tt),yVals,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        delete(currentTTFHandle)
        currentTTFHandle = semilogx(freqDomain,model(freqDomain,psiParamsQuest),'-r');
        legend('Veridical model','Stimulus Outcomes','Best Fit from Q+','Location','northwest');
        
        % Entropy by trial
        subplot(3,1,3)
        delete(currentEntropyHandle)
        entropyAfterTrial(1:tt)=questData.entropyAfterTrial;
        plot(1:nTrials,entropyAfterTrial,'*k');
        xlim([1 nTrials]);
        ylim([0 nanmax(entropyAfterTrial)]);
        xlabel('Trial number');
        ylabel('Entropy');
        drawnow
        
    end
    
    
end

%% Adjust Beta/MaxBOLD tradeoff
% For the final parameter estimate, we want to assume that Beta is 1 and
% maxBOLD is whatever it WOULD BE if beta were 1. 

% Grab our current beta estimate is: 
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
betaGuess = psiParamsQuest(betaIndex);

% Divide maxBOLD by our beta estimate: (beta / beta) = 1, so
% new maxBOLD = maxBOLD/beta 
maxBOLD = maxBOLD.*betaGuess;

% Now run through the fitting steps again with the new maxBOLD
thePacket = createPacket('nTrials',tt,...,
    'trialLengthSecs',trialLength,...,
    'stimulusStructDeltaT',stimulusStructDeltaT);

[outcomes] = ...
    tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
    'maxBOLDSimulated',maxBOLDSimulated,...
    'rngSeed',rngSeed.Seed,...,
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
    'lowerBounds', lowerBounds,'upperBounds',upperBounds,...
    'plausibleLowerBounds',lowerBounds,'plausibleUpperBounds',upperBounds);


fprintf('\nMax BADS parameters:               '); 
for i = 1:length(psiParamsFit)
    fprintf('%s: %0.3f ',paramNamesInOrder{i},psiParamsFit(i));
end

fprintf('\nmaxBOLD estimate: %0.3f\n',maxBOLD);

if showPlots
    figure('Position', [10 10 figWidth figHeight]);
    hold on;
    semilogx(freqDomain,model(freqDomain,simulatedPsiParams),'-k');
    semilogx(freqDomain,model(freqDomain,psiParamsQuest),'-r');
    semilogx(freqDomain,model(freqDomain,psiParamsFit),'-b');
    set(gca,'XScale', 'log')
    xlabel('log stimulus Frequency [Hz]');
    ylabel('Relative response amplitude');
    legend('Veridical model','Best Fit from Q+','Best Fit from BADS','Location','northwest');
end

%% Output
T = array2table(psiParamsFit,'VariableNames',paramNamesInOrder);
T.maxBOLD = maxBOLD;
outfilename = horzcat(outNum,'.csv');
%save(outfilename,psiParamsFit);
writetable(T,outfilename);

end




