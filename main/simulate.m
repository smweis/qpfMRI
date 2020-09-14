function [qpfmriResults]=simulate(myQpfmriParams, myQpParams, varargin)
%% [qpfmriResults]=simulate(myQpfmriParams, varargin)
% A script that will simulate fMRI BOLD data and fit a model with or
% without Q+ control
%
% Syntax:
%  [qpfmriResults]=simulate(myQpfmriParams, varargin)
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
% Optional key/value pairs (used in plotting):
%   'showPlots'             - Logical: (Default = false)
%                             Whether to show plots.
%   'saveFigs'              - Logical (Default = false)
%                             Whether to save the figures as output.
%   'saveGif'               - Logical (Default = false)
%                             Whether to save the animated plot as a gif.
% Outputs:
%   qpfmriResults         - Struct. Results of simulation. 
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

% Set the number of outcome categories / bins.
nOutcomes = 15;

% Do you want to see plots?
showPlots = true; 

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


% Run the simulation. 
[qpfmriResults]=simulate(myQpfmriParams,myQpParams,'showPlots',showPlots);

%}

%% TO DO
%1. Allow modification of baseline and maxBOLD trials

%% Parse Inputs
p = inputParser;

p.addRequired('myQpfmriParams',@isstruct);
p.addRequired('myQpParams',@isstruct);

% Optional params for plotting
p.addParameter('showPlots',false,@islogical);
p.addParameter('saveFigs',false,@islogical);
% Parse inputs
p.parse( myQpfmriParams, myQpParams, varargin{:});


%% Initialize struct to save results out.
qpfmriResults = myQpfmriParams;

% Initialize maxBOLD
qpfmriResults.maxBOLDoverTrials = nan(1,myQpfmriParams.nTrials);
maxBOLDLatestGuess = myQpfmriParams.maxBOLDInitialGuess;
qpfmriResults.psiParamsQuest = nan(myQpfmriParams.nTrials,length(fieldnames(myQpfmriParams.paramsDomain)));
qpfmriResults.entropyOverTrials = cell(1,myQpfmriParams.nTrials);

% Set up save info and directory
folderName = ['.' filesep myQpfmriParams.outFolder];
fileName = ['sim_' myQpfmriParams.outNum '.mat'];
if ~exist(folderName,'dir')
    mkdir(folderName);
end


%% Handle seed

% If seed is a keyword, pick a random seed. We use shuffle so it's new
% every time.
if strcmp(myQpfmriParams.seed,'choose')
    rngSeed = rng('shuffle'); rngSeed = rng('shuffle');
else 
    fprintf('Initial random call %.04f\n',rand);
    rngSeed = rng(myQpfmriParams.seed); rngSeed = rng(myQpfmriParams.seed);
    fprintf('Random call after seed %.04f\n',rand);
end



%% Derive some lower and upper bounds from the parameter ranges. This is
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
% Warn the user that initializing has started and may take a minute.
tic
fprintf('Initializing Q+. This may take a minute...\n');
questData = qpInitialize(myQpParams);
toc

% Tack on a continuous output simulated observer to myQpParams
myQpParams.continuousPF = @(f) myQpfmriParams.model(f,myQpfmriParams.simulatedPsiParams);

% Calculate the lower headroom bin offset. We'll use this for plotting.
nLower = max([1 round(myQpfmriParams.headroom*myQpParams.nOutcomes)]);
nUpper = max([1 round(myQpfmriParams.headroom*myQpParams.nOutcomes)]);
nMid = myQpParams.nOutcomes - nLower - nUpper;

% Create a copy of Q+
questDataUntrained = questData;

% Create a stimulusVec to hold the trial across the loops
qpfmriResults.stimulusVec = nan(1,myQpfmriParams.nTrials);
% This allows us to store how each trial was generated
qpfmriResults.stimulusVecTrialTypes = cell(1,myQpfmriParams.nTrials);

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

    % First trials require baseline or maxBOLD. 
    if tt<=6
        % Baseline trial
        if mod(tt,2) > 0 
            qpfmriResults.stimulusVec(tt) = myQpfmriParams.baselineStimulus;
            qpfmriResults.stimulusVecTrialTypes{tt} = 'baseline';
        % maxBOLD trial
        else 
            qpfmriResults.stimulusVec(tt) = myQpfmriParams.maxBOLDStimulus;
            qpfmriResults.stimulusVecTrialTypes{tt} = 'maxBOLD';
        end
    % Optionally, every X trials can alternate as baseline or maxBOLD. 
    % Baseline trial
    elseif mod(tt,5) == 0
        qpfmriResults.stimulusVec(tt) = myQpfmriParams.baselineStimulus;
        qpfmriResults.stimulusVecTrialTypes{tt} = 'baseline';
    % maxBOLD trial
    elseif mod(tt,10) == 0
        qpfmriResults.stimulusVec(tt) = myQpfmriParams.maxBOLDStimulus;
        qpfmriResults.stimulusVecTrialTypes{tt} = 'maxBOLD';
    % Every other trial will be selected randomly, or by Q+
    else
        if ~myQpfmriParams.qpPres
            qpfmriResults.stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            qpfmriResults.stimulusVecTrialTypes{tt} = 'random';
        else
            qpfmriResults.stimulusVec(tt) = qpQuery(questData);
            qpfmriResults.stimulusVecTrialTypes{tt} = 'qplus';
        end
    end
    

    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if tt > 2
        maxBOLDLatestGuess = maxBOLDLatestGuess.*qpfmriResults.psiParamsQuest(tt-1,myQpfmriParams.betaIndex);
    end

    % Create a packet
    thePacket = createPacket(myQpfmriParams,tt);

    % Obtain outcomes from tfeUpdate 
    [outcomes, modelResponseStruct, thePacketOut, ~, baselineEstimate] = ...
        tfeUpdate(thePacket, myQpParams, myQpfmriParams, qpfmriResults.stimulusVec, maxBOLDLatestGuess, 'rngSeed',rngSeed);

    % Grab a naive copy of questData
    questData = questDataUntrained;

    % Update quest data structure. This is the slow step in the simulation.
    for yy = 1:tt
        questData = qpUpdate(questData,qpfmriResults.stimulusVec(yy),outcomes(yy));
    end
    
    % Save individual run output.
    psiParamsIndex = qpListMaxArg(questData.posterior);
    qpfmriResults.psiParamsQuest(tt,:) = questData.psiParamsDomain(psiParamsIndex,:);
    qpfmriResults.maxBOLDoverTrials(tt) = maxBOLDLatestGuess;
    
    %% A QUESTION HERE! Do we want to store the MIN entropy for each trial, 
    % or do we want to store ALL entropy values (1 per trial) for each trial
    qpfmriResults.entropyOverTrials{tt} = questData.entropyAfterTrial;
    
    %% 
    % Calculate BOLD values from outcomes, scaled to maxBOLD and the
    % baseline estimate.
    yVals = (outcomes - nLower - 1)./nMid;
    yVals = yVals*qpfmriResults.psiParamsQuest(tt,myQpfmriParams.betaIndex)*maxBOLDLatestGuess;
    yValsPlusBaseline = yVals + baselineEstimate;

    % Print results to the console.
    trialString = sprintf('\nTrial %d: %s Stimulus Selection\n',tt,qpfmriResults.stimulusVecTrialTypes{tt});
    fprintf(trialString);
    stimString = sprintf('Stimulus: %0.3f\n',qpfmriResults.stimulusVec(tt));
    fprintf(stimString);
    resultString = sprintf('Output value %.03f\n',yVals(end));
    fprintf(resultString);
    fprintf('\nQ+ parameters\n');
    for i = 1:length(qpfmriResults.psiParamsQuest(tt,:))
        fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},qpfmriResults.psiParamsQuest(tt,i));
    end
    fprintf('\nmaxBOLD estimate = %0.3f\n',maxBOLDLatestGuess);
    fprintf('\n');
    
    %% Plot the ongoing results
    % Update the plots
    if p.Results.showPlots
        % Draw plots
        [mainFig,handleStruct] = drawPlots(myQpfmriParams,myQpParams,...,
            qpfmriResults,mainFig,handleStruct,thePacketOut,modelResponseStruct,...,
            yVals,yValsPlusBaseline);
    end
    
end

%% Adjust Beta/MaxBOLD tradeoff
% For the final parameter estimate, we want to assume that Beta is 1 and
% maxBOLD is whatever it WOULD BE if beta were 1. 

% Grab our current beta estimate: 
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,qpfmriResults.psiParamsQuest(tt,:),questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds,...
    'plausibleLowerBounds',lowerBounds,'plausibleUpperBounds',upperBounds);
maxBOLDLatestGuess = maxBOLDLatestGuess.*psiParamsFit(myQpfmriParams.betaIndex);

% Now run through the fitting steps again with the new maxBOLD
thePacket = createPacket(myQpfmriParams,myQpfmriParams.nTrials);
[outcomes, ~, ~, ~, ~] = ...
        tfeUpdate(thePacket, myQpParams, myQpfmriParams, qpfmriResults.stimulusVec, maxBOLDLatestGuess, 'rngSeed',rngSeed);
questData = questDataUntrained;
for yy = 1:tt
    questData = qpUpdate(questData,qpfmriResults.stimulusVec(yy),outcomes(yy));
end

%% Print some final output to the log
fprintf('--------------------------------------------------------------\n');
fprintf('FINAL VALUES\n');
fprintf('--------------------------------------------------------------\n');
% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questData.posterior);
qpfmriResults.psiParamsQuestFinal = questData.psiParamsDomain(psiParamsIndex,:);
% Print simulated parameters.
fprintf('Simulated parameters:              '); 
for i = 1:length(myQpfmriParams.simulatedPsiParams)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},myQpfmriParams.simulatedPsiParams(i));
end
% Print Q+ Fit parameters.
fprintf('\nMax posterior QUEST+ parameters:   '); 
for i = 1:length(qpfmriResults.psiParamsQuestFinal)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},qpfmriResults.psiParamsQuestFinal(i));
end

% Prepare for BADS fit.
psiParamsBads = qpfmriResults.psiParamsQuestFinal;
psiParamsBads(myQpfmriParams.betaIndex) = 1;

% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
qpfmriResults.psiParamsBadsFinal = qpFitBads(questData.trialData,questData.qpPF,psiParamsBads,questData.nOutcomes,...
    'lowerBounds', lowerBoundsConstrained,'upperBounds',upperBoundsConstrained,...
    'plausibleLowerBounds',lowerBoundsConstrained,'plausibleUpperBounds',upperBoundsConstrained);
% Print BADS best fit parameters.
fprintf('\nMax BADS parameters:               '); 
for i = 1:length(qpfmriResults.psiParamsBadsFinal)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},qpfmriResults.psiParamsBadsFinal(i));
end

% Print maxBOLD
fprintf('\nmaxBOLD estimate: %0.3f\n',maxBOLDLatestGuess);

qpfmriResults.maxBOLDFinal = maxBOLDLatestGuess;

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



save(fullfile(folderName,fileName),'qpfmriResults');


end




