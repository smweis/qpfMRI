function [qpfmriResults,questDataCopy]=simulate(myQpfmriParams, myQpParams, varargin)
%% [qpfmriResults,questDataCopy]=simulate(myQpfmriParams, varargin)
% A script that will simulate fMRI BOLD data and fit a model with or
% without Q+ control
%
% Syntax:
%  [qpfmriResults,questDataCopy]=simulate(myQpfmriParams, varargin)
%
% Description:
%	Takes in a model and a possible set of parameters and whether or not Q+ 
%   is in control of things flag. 
%
% Inputs:
% Required Inputs
%   myQpfmriParams          - Struct. Set of parameters used for qpfmri.
%                             See qpfmriParams function for more details
%   myQpParams              - Struct. Set of parameters used for qp.
%                             See qpParams function for more details
% Optional key/value pairs:
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
%   'saveFigs'              - Logical (Default = false)
%                             Whether to save the figures as output.
%   'saveGif'               - Logical (Default = false)
%                             Whether to save the animated plot as a gif.
% Outputs:
%   qpfmriResults         - Struct. Results of simulation. 
%   questDataCopy         - Struct. Copy of initialized questData.
%
%Examples: 
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
nTrials = 30;

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

[myQpfmriParams,myQpParams] = qpfmriParams(model,paramsDomain,'qpPres',qpPres,...,
'stimulusDomain',stimulusDomain,'stimulusDomainSpacing',stimulusDomainSpacing,...,
'noiseSD',noiseSD,'nTrials',nTrials,'maxBOLDSimulated',maxBOLDSimulated,...,
'trialLength',trialLength,'nOutcomes',nOutcomes);


% Note, this will save a copy of questData after it is initialized. 
[qpfmriResults,questDataCopy]=simulate(myQpfmriParams,myQpParams,'showPlots',showPlots);
---------------------------------------------------------------------------
Example 1 Time Saver code:

% Time saver for debugging: After running one of the above examples, keep
% everything in memory and run the line below. Especially useful if the
% paramsDomain is large or multi-dimensional.

% If paramsDomain is not changed, the following line can be run with
% questDataCopy as an optional argument to save the initialization step.
[qpfmriResults,questDataCopy]=simulate(myQpfmriParams,'showPlots',showPlots,...,
'questDataCopy',questDataCopy);

%}

%% TO DO
%1. Separate out plotting as separate functions. 
%2. Create results struct to save out. 


%% Parse Inputs
p = inputParser;

p.addRequired('myQpfmriParams',@isstruct);
p.addRequired('myQpParams',@isstruct);
% Check if there's a copy of questData already
p.addParameter('questDataCopy',{},@isstruct);

% Optional params for plotting
p.addParameter('showPlots',false,@islogical);
p.addParameter('saveFigs',false,@islogical);
% Parse inputs
p.parse( myQpfmriParams, myQpParams, varargin{:});


% Initialize struct to save results out.
qpfmriResults = myQpfmriParams;

% Initialize maxBOLD
qpfmriResults.maxBOLDoverTrials = nan(1,myQpfmriParams.nTrials);
maxBOLDLatestGuess = myQpfmriParams.maxBOLDInitialGuess;

%% Handle inputs

% If seed is a keyword, pick a random seed. We use shuffle so it's new
% every time.
if strcmp(myQpfmriParams.seed,'choose')
    rngSeed = rng('shuffle'); rngSeed = rng('shuffle');
else 
    fprintf('Initial random call %.04f\n',rand);
    rngSeed = rng(myQpfmriParams.seed); rngSeed = rng(myQpfmriParams.seed);
    fprintf('Random call after seed %.04f\n',rand);
end


%% Select the veridical psychometric paramteers for the function.
% Pick some random params to simulate if none provided but set beta to 1.
% We require the simulated parameters to result in a baseline trial = 0.
if isempty(myQpfmriParams.simulatedPsiParams)
    myQpfmriParams.simulatedPsiParams = zeros(1,length(myQpfmriParams.paramNamesInOrder));
    stillSearching = true;
    while stillSearching
        for i = 1:length(myQpfmriParams.paramNamesInOrder)
            myQpfmriParams.simulatedPsiParams(i) = randsample(myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{i}),1);
        end
        
        % Beta is always one for simulations
        myQpfmriParams.simulatedPsiParams(myQpfmriParams.betaIndex) = 1;
        
        % Simulated noise is selected from a random sample of noiseSD
        myQpfmriParams.simulatedPsiParams(myQpfmriParams.sigmaIndex) = randsample(myQpfmriParams.noiseSD,1);
        
        if abs(myQpfmriParams.model(myQpfmriParams.baselineStimulus,myQpfmriParams.simulatedPsiParams)) < myQpfmriParams.simulatedPsiParams(myQpfmriParams.betaIndex)/10000 && ...
                abs(myQpfmriParams.model(myQpfmriParams.maxBOLDStimulus,myQpfmriParams.simulatedPsiParams)) < 1 && ...
                abs(myQpfmriParams.model(myQpfmriParams.maxBOLDStimulus,myQpfmriParams.simulatedPsiParams)) > .99
            stillSearching = false;
        end
    end
    

else
    myQpfmriParams.simulatedPsiParams = zeros(1,length(myQpfmriParams.paramNamesInOrder));
    
    % A relic of naming simulatedPsiParams noise parameter as sigma.
    for i = 1:length(myQpfmriParams.paramNamesInOrder)
        if i == myQpfmriParams.sigmaIndex && ~isfield(myQpfmriParams.simulatedPsiParams,'sigma')
            myQpfmriParams.simulatedPsiParams(myQpfmriParams.sigmaIndex) = noiseSD;
        else
            myQpfmriParams.simulatedPsiParams(i) =  myQpfmriParams.simulatedPsiParams.(myQpfmriParams.paramNamesInOrder{i});
        end
    end
    
    % Beta will converge to 1 as maxBOLD gets closer and closer to the
    % simulated maxBOLD. As a result, when simulating data, beta should always
    % be set to 1. 
    myQpfmriParams.simulatedPsiParams(myQpfmriParams.betaIndex) = 1;
    assert(myQpfmriParams.simulatedPsiParams(myQpfmriParams.betaIndex)==1,'Simulated Beta should always be 1.');
    
    % Select simulated psychometric parameters whose range of outputs
    % are between 0 and 1 for the model being used. 
    if abs(myQpfmriParams.model(myQpfmriParams.baselineStimulus,myQpfmriParams.simulatedPsiParams)) < myQpfmriParams.simulatedPsiParams(myQpfmriParams.betaIndex)/10000
        warning('Simulated psychometric parameters will result in minimum values below 0.\nMin possible value = %.02f',abs(myQpfmriParams.model(myQpfmriParams.baselineStimulus,myQpfmriParams.simulatedPsiParams)));
    elseif abs(myQpfmriParams.model(myQpfmriParams.baselineStimulus,myQpfmriParams.simulatedPsiParams)) > .01
        warning('Simulated psychometric parameters will result in minimum values greater than 0.\nMin possible value = %.02f',abs(myQpfmriParams.model(myQpfmriParams.baselineStimulus,myQpfmriParams.simulatedPsiParams)));
    elseif abs(myQpfmriParams.model(myQpfmriParams.maxBOLDStimulus,myQpfmriParams.simulatedPsiParams)) < .99
        warning('Simulated psychometric parameters will result in maximum BOLD values below 1.\nMax possible value = %.02f',myQpfmriParams.model(myQpfmriParams.maxBOLDStimulus,myQpfmriParams.simulatedPsiParams));
    end
end

% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = zeros(1,length(myQpfmriParams.paramNamesInOrder));
upperBounds = zeros(1,length(myQpfmriParams.paramNamesInOrder));
for i = 1:length(myQpfmriParams.paramNamesInOrder)
    lowerBounds(i) = myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{i})(1);
    upperBounds(i) = myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{i})(end);
end

% Constrain bounds on beta to be very tight around 1.
lowerBoundsConstrained = lowerBounds;
upperBoundsConstrained = upperBounds;
lowerBoundsConstrained(myQpfmriParams.betaIndex) = .999;
upperBoundsConstrained(myQpfmriParams.betaIndex) = 1.001;

% Make sure 1 is a member of the beta parameter domains. If it's not, add
% it in. 
if ~ismember(1,myQpfmriParams.paramsDomain.beta)
    warning('The domain for beta should always include 1. Adding 1 to beta parameter domain. This may create non-linear spacing.');
    myQpfmriParams.paramsDomain.beta = sort(myQpfmriParams.paramsDomain.beta);
    myQpfmriParams.paramsDomain.beta = [myQpfmriParams.paramsDomain.beta(myQpfmriParams.paramsDomain.beta<1) 1 myQpfmriParams.paramsDomain.beta(myQpfmriParams.paramsDomain.beta>1)];
end

% Veridical values should be within domain bounds, except for sigma.
for param = 1:length(myQpfmriParams.paramNamesInOrder)
    if ~strcmp(myQpfmriParams.paramNamesInOrder{param},'sigma')
        assert(lowerBounds(param) <= myQpfmriParams.simulatedPsiParams(param)...
            && upperBounds(param) >= myQpfmriParams.simulatedPsiParams(param),...,
            'Parameter %s is not within the bounds of the parameter domain.',myQpfmriParams.paramNamesInOrder{param});
    end
end

% Create a simulated observer with binned output
myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,myQpfmriParams.simulatedPsiParams);

%% Initialize Q+. 
% Save some time if Q+ has already been initialized with the same
% parameters and the user has passed it as an input. 

% Time saver will use a blank copy of questData if it's been initialized
% and passed in as an argument. 
try
    questDataCopy = p.Results.questDataCopy;
catch
    warning('No questData found. Will initialize, which will take some time.');
end

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
nLower = max([1 round(myQpfmriParams.headroom*myQpParams.nOutcomes)]);
nUpper = max([1 round(myQpfmriParams.headroom*myQpParams.nOutcomes)]);
nMid = myQpParams.nOutcomes - nLower - nUpper;

% Tack on a continuous output simulated observer to myQpParams
myQpParams.continuousPF = @(f) myQpfmriParams.model(f,myQpfmriParams.simulatedPsiParams);

% Create a copy of Q+
questDataUntrained = questData;

% Create a stimulusVec to hold the trial across the loops
stimulusVec = nan(1,myQpfmriParams.nTrials);
% This allows us to store how each trial was generated
stimulusVecTrialTypes = cell(1,myQpfmriParams.nTrials);
%% Plotting features

% Create a plot in which we can track the model progress
if p.Results.showPlots
    [mainFig,handleStruct] = initializePlots(myQpfmriParams,myQpParams);
end

%% Print the simulated psychometric parameters:
fprintf('Simulated parameters:\n'); 
for i = 1:length(myQpfmriParams.simulatedPsiParams)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},myQpfmriParams.simulatedPsiParams(i));
end
fprintf('\n'); 

%% Run simulated trials
for tt = 1:myQpfmriParams.nTrials

    % If it is the first three trials we force a baseline or maxBOLD event
    % Require 3 each. 
    if tt<=6
        if mod(tt,2) > 0 % Baseline trial
            stimulusVec(tt) = myQpfmriParams.baselineStimulus;
            stimulusVecTrialTypes{tt} = 'baseline';
        else % maxBOLD trial
            stimulusVec(tt) = myQpfmriParams.maxBOLDStimulus;
            stimulusVecTrialTypes{tt} = 'maxBOLD';
        end
    % After that, every 5th trial will be a baseline or a maxBOLD trial.    
    elseif mod(tt,5) == 0 % Baseline trial
        stimulusVec(tt) = myQpfmriParams.baselineStimulus;
        stimulusVecTrialTypes{tt} = 'baseline';
    elseif mod(tt,10) == 0 % maxBOLD trial
        stimulusVec(tt) = myQpfmriParams.maxBOLDStimulus;
        stimulusVecTrialTypes{tt} = 'maxBOLD';
    % Every other trial will be selected randomly, or by Q+
    else
        if ~myQpfmriParams.qpPres
            % get random stimulus
            stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            stimulusVecTrialTypes{tt} = 'random';
        else
            % get next stimulus from Q+ 
            stimulusVec(tt) = qpQuery(questData);
            stimulusVecTrialTypes{tt} = 'qplus';
        end
    end
    
    trialString = sprintf('\nTrial %d: %s Stimulus Selection\n',tt,stimulusVecTrialTypes{tt});
    fprintf(trialString);
    stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(tt));
    fprintf(stimString);
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % that could be evoked by a stimulus (relative to the baseline
    % stimulus), which is the beta value of the model
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if tt > 2
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        % Update maxBOLD by the current best guess for Beta.
        maxBOLDLatestGuess = maxBOLDLatestGuess.*psiParamsQuest(myQpfmriParams.betaIndex);
    end

    % Create a packet
    thePacket = createPacket(myQpfmriParams,tt);

    % Obtain outcomes from tfeUpdate 
    [outcomes, modelResponseStruct, thePacketOut, ~, baselineEstimate] = ...
        tfeUpdate(thePacket, myQpParams, myQpfmriParams, stimulusVec, maxBOLDLatestGuess, 'rngSeed',rngSeed);

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
    yVals = yVals*psiParamsQuest(myQpfmriParams.betaIndex)*maxBOLDLatestGuess;
    yValsPlusBaseline = yVals + baselineEstimate;
    
    % Print results to the console.
    resultString = sprintf('Output value %.03f\n',yVals(end));
    fprintf(resultString);
    fprintf('\nQ+ parameters\n');
    for i = 1:length(psiParamsQuest)
        fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},psiParamsQuest(i));
    end
    fprintf('\nmaxBOLD estimate = %0.3f\n',maxBOLDLatestGuess);
    fprintf('\n');
    
    %% Plot the ongoing results
    % Update the plots
    if p.Results.showPlots

        % Draw plots
        [mainFig,handleStruct] = drawPlots(myQpfmriParams,myQpParams,...,
            stimulusVec,mainFig,handleStruct,thePacketOut,modelResponseStruct,...,
            yVals,yValsPlusBaseline,psiParamsQuest,maxBOLDLatestGuess,questData.entropyAfterTrial,'trialTypes',stimulusVecTrialTypes);
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
maxBOLDLatestGuess = maxBOLDLatestGuess.*psiParamsFit(myQpfmriParams.betaIndex);

% Now run through the fitting steps again with the new maxBOLD
thePacket = createPacket(myQpfmriParams,myQpfmriParams.nTrials);

[outcomes, ~, ~, ~, ~] = ...
        tfeUpdate(thePacket, myQpParams, myQpfmriParams, stimulusVec, maxBOLDLatestGuess, 'rngSeed',rngSeed);

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
for i = 1:length(myQpfmriParams.simulatedPsiParams)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},myQpfmriParams.simulatedPsiParams(i));
end
% Print Q+ Fit parameters.
fprintf('\nMax posterior QUEST+ parameters:   '); 
for i = 1:length(psiParamsQuest)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},psiParamsQuest(i));
end
% Prepare for BADS fit.
psiParamsBads = psiParamsQuest;
psiParamsBads(myQpfmriParams.betaIndex) = 1;

% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsBads,questData.nOutcomes,...
    'lowerBounds', lowerBoundsConstrained,'upperBounds',upperBoundsConstrained,...
    'plausibleLowerBounds',lowerBoundsConstrained,'plausibleUpperBounds',upperBoundsConstrained);
% Print BADS best fit parameters.
fprintf('\nMax BADS parameters:               '); 
for i = 1:length(psiParamsFit)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},psiParamsFit(i));
end
% Print maxBOLD
fprintf('\nmaxBOLD estimate: %0.3f\n',maxBOLDLatestGuess);


%% Output
% Save plots as pdf

% Save figures
if p.Results.showPlots
    if p.Results.saveFigs
        delete(findall(gcf,'type','annotation'));
        set(mainFig,'Units','Inches');
        pos = get(mainFig,'Position');
        set(mainFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        print(mainFig,'simulate.pdf','-dpdf','-r0');
    
        set(paramsFig,'Units','Inches');
        pos = get(paramsFig,'Position');
        set(paramsFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        print(paramsFig,'paramsFigs.pdf','-dpdf','-r0');
    
        set(finalFig,'Units','Inches');
        pos = get(finalFig,'Position');
        set(finalFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        print(finalFig,'finalFig.pdf','-dpdf','-r0');

    end
end

% Save results in struct
resultsOut = array2table(psiParamsFit,'VariableNames',paramNamesInOrder);
resultsOut.maxBOLD = maxBOLD;
resultsFileName = [outNum 'Results.csv'];
resultsFolderName = ['.' filesep myQpfmriParams.outFolder filesep 'results'];

if ~exist(resultsFolderName,'dir')
    mkdir(resultsFolderName);
end


resultsOut = array2table(psiParamsFit,'VariableNames',paramNamesInOrder);
resultsOut.maxBOLDLatestGuess = maxBOLDLatestGuess;

writetable(resultsOut,[resultsFolderName filesep resultsFileName]);

simulatedParamsOut = array2table(myQpfmriParams.simulatedPsiParams,'VariableNames',paramNamesInOrder);
simulatedParamsOut.maxBOLDSimulated = myQpfmriParams.maxBOLDSimulated;
simulatedParamsOut.nOutcomes = myQpfmriParams.nOutcomes;

paramsFileName = [outNum 'Params.csv'];
paramsFolderName = ['.' filesep myQpfmriParams.outFolder filesep 'params'];

if ~exist(paramsFolderName,'dir')
    mkdir(paramsFolderName);
end

writetable(simulatedParamsOut,[paramsFolderName filesep paramsFileName]);

end




