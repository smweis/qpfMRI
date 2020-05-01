function [modelResponseStruct,thePacketOut,questDataCopy]=validate_qpDoETFE_simulate(model_params, control_params, sim_type, varargin)
%% A validation script to evaluate how the DoE model is working
%
% Syntax:
%  [modelResponseStruct,thePacketOut,questDataCopy]=validate_qpDoETFE_simulate(model_params, control_params, sim_type, debug)
%
% Description:
%	Takes in simulated parameters for the model, the control (TR, trial
%	length), and whether or not Q+ is in control of things, plus a debug
%	flag. 
%
% Inputs:
%   model_params          - 1x5 vector. sr, kr1, kr2, beta, sigma
%
%   control_params        - 1x2 vector. TR (secs), trial length (msecs)
%
%   sim_type              - Logical. 1 for Q+, 0 for random.
%
% Optional key/value pairs (used in fitting):
%  'questDataCopy'        - Struct. Q+ data structure.
%
% Outputs:
%   modelResponseStruct   - Struct. The simulated response struct from tfe 
%                           method fitResponse
%   thePacketOut          - Struct. The updated packet with response struct
%                           completed.
%   questDataCopy         - Struct. Q+ data structure.
%
%
% Examples:
%{

% SIMULATED BETA MUST ALWAYS BE UNITY
model_params = [1.05 .01 .06 1.00 .4]; 
control_params = [800 12]; %TR (secs), trial length (msecs)
sim_type = false; %Q+ (if true), random (if false)
[modelResponseStruct,thePacketOut,questDataCopy]=validate_qpDoETFE_simulate(model_params, control_params, sim_type);


% To run again and skip reinitializing (only do so if you have NOT changed the parameter domains):

[modelResponseStruct,thePacketOut,questDataCopy]=validate_qpDoETFE_simulate(model_params, control_params, sim_type,'questDataCopy',questDataCopy);
%}



%% Are we debugging?
p = inputParser;

% Required input
p.addRequired('model_params',@isvector);
p.addRequired('control_params',@isvector);
p.addRequired('sim_type',@islogical);


% Optional params used in fitting
p.addParameter('questDataCopy', @isstruct);

% Parse
p.parse(model_params, control_params, sim_type, varargin{:});

try
    questDataCopy = p.Results.questDataCopy;
catch
    warning('No questData found. Will initialize, which will take some time.');
end

%% Are we simulating old fashioned constant stimuli?
simulateConstantStimuli = sim_type; 


%% Define the veridical model params

% Leave the simulatedPsiParams empty to try a random set of params.
% Here are some params to try for specific TTF shapes:
%  A low-pass TTF in noisy fMRI data: [10 1 0.83 1]
%  A band-pass TTF in noisy fMRI data: [1.47 1.75 0.83 1]
%simulatedPsiParams = [4 1 1 1 0];
%simulatedPsiParams = [0.9998 0.0132 0.7755 1 0];
simulatedPsiParams = model_params;

% Some information about the trials?
nTrials = 30; % how many trials
trialLengthSecs = control_params(2); % seconds per trial (12)
stimulusStructDeltaT = 100; % the resolution of the stimulus struct in msecs

% True size of the BOLD response
maxBOLDSimulated = 1.5;

% Initial guess for the max size of the evoked BOLD response
maxBOLD = 1.0;

% Which stimulus (in freq Hz) is the "baseline" stimulus? This stimulus
% should be selected with the expectation that the neural response to this
% stimulus will be minimal as compared to all other stimuli.
baselineStimulus = 0;

% How talkative is the simulation?
showPlots = true;
verbose = true;


%% Set up Q+

% Get the default Q+ params
myQpParams = qpParams;

% Add the stimulus domain. ~Log spaced frequencies between 0 and 30 Hz
myQpParams.stimParamsDomainList = {[baselineStimulus,1.875,3.75,7.5,15,30]};
nStims = length(myQpParams.stimParamsDomainList{1});

% The number of outcome categories.
myQpParams.nOutcomes = 51;

% The headroom is the proportion of outcomes that are reserved above and
% below the min and max output of the DoE model to account for noise
headroom = 0.1;

% Create an anonymous function from qpDoETemporalModel in which we
% specify the number of outcomes for the y-axis response
myQpParams.qpPF = @(f,p) qpDoETemporalModel(f,p,myQpParams.nOutcomes,headroom);

% Define the parameter ranges
% NOTE: IF YOU CHANGE THESE, YOU MUST RE-INITIALIZE Q+
Sr = 0.899:0.025:1.099;
k1 = 0.001:0.0005:0.01;
k2 = 0.001:0.01:.2;
beta = 0.4:0.2:2; % Amplitude of the scaled response; should converge to unity
    % NEED TO ENSURE THAT THE VALUE OF UNITY IS PRESENT IN THE SET OF
    % BETA VALUES
sigma = 0.3:0.2:1;	% Standard deviation of the scaled (0-1) noise
myQpParams.psiParamsDomainList = {Sr, k1, k2, beta, sigma};

% Pick some random params to simulate if none provided (but set the neural
% noise to zero)
if isempty(simulatedPsiParams)
    simulatedPsiParams = [randsample(Sr,1) randsample(k1,1) randsample(k2,1) randsample(beta,1) 0];
end

% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = [Sr(1) k1(1) k2(1) beta(1) sigma(1)];
upperBounds = [Sr(end) k1(end) k2(end) beta(end) sigma(end)];

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
myQpParams.continuousPF = @(f) doeTemporalModel(f,simulatedPsiParams);


% Create a full length packet
thePacket = createPacket('nTrials',nTrials,...,
    'trialLengthSecs',trialLengthSecs,...,
    'stimulusStructDeltaT',stimulusStructDeltaT);
 

% Create a plot in which we can track the model progress
if showPlots
    figure
    
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
    freqDomain = logspace(log10(0.01),log10(100),100);
    semilogx(freqDomain,doeTemporalModel(freqDomain,simulatedPsiParams),'-k');
    ylim([-0.5 1.5]);
    xlabel('log stimulus Frequency [Hz]');
    ylabel('Relative response amplitude');
    title('Estimate of DoE TTF');
    hold on
    currentOutcomesHandle = scatter(nan,nan);
    currentTTFHandle = plot(freqDomain,doeTemporalModel(freqDomain,simulatedPsiParams),'-k');
    
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


% Create and save an rng seed to use for this simulation
rngSeed = rng();

% Create a copy of Q+
questDataUntrained = questData;

% Create a stimulusVec to hold the trial across the loops
stimulusVec = nan(1,nTrials);

%% Run simulated trials
for tt = 1:nTrials
    
    % Ask QP to supply our next stimulus. If it is the first two trials
    % we force a baseline event
    if tt<=2
        stimulusVec(tt) = baselineStimulus;
    else
        if simulateConstantStimuli
            % get random stimulus
            stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
        else
            % get next stimulus from Q+
            stimulusVec(tt) = qpQuery(questData);
        end
    end
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % that could be evoked by a stimulus (relative to the baseline
    % stimulus), which is the beta value of the model
    try 
        psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds)
        maxBOLD = maxBOLD.*psiParamsFit(4);
    catch
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        maxBOLD = maxBOLD.*psiParamsQuest(4);
    end

    % Create a packet
    thePacket = createPacket('nTrials',tt,...,
        'trialLengthSecs',trialLengthSecs,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);
    
    % If we were not in simulation mode, we would add the BOLD fMRI
    % response to the packet here. Instead, it will be simulatd within
    % tfeUpdate.
    
    % Obtain outcomes from tfeUpdate 
    [outcomes, modelResponseStruct, thePacketOut] = ...
        tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
        'maxBOLDSimulated',maxBOLDSimulated,...
        'rngSeed',rngSeed,...,
        'maxBOLD',maxBOLD,...,
        'TRmsecs', control_params(1));
   
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
        stimulusVecPlot(stimulusVecPlot==0)=0.01;
        delete(currentOutcomesHandle);
        currentOutcomesHandle = scatter(stimulusVecPlot(1:tt),yVals,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        delete(currentTTFHandle)
        currentTTFHandle = semilogx(freqDomain,doeTemporalModel(freqDomain,psiParamsQuest),'-r');
        legend('Veridical','Stimulus Outcomes','Best Fit from Q+');
    
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

% Done with the simulation
if verbose
    fprintf('\n');
end



% Remove some of the larger variables for output
%questDataOut = questData;
%questDataOut.precomputedOutcomeProportions = [];
%questDataOut.psiParamsDomain = [];

%% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
fprintf('Simulated parameters:              %0.3f, %0.3f, %0.3f, %0.3f, %0.3f \n', ...
    simulatedPsiParams(1),simulatedPsiParams(2),simulatedPsiParams(3),simulatedPsiParams(4),simulatedPsiParams(5));
fprintf('Max posterior QUEST+ parameters:   %0.3f, %0.3f, %0.3f, %0.3f, %0.3f \n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));

%% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds,'diagnostics','on','display','on');
fprintf('Maximum likelihood fit parameters: %0.3f, %0.3f, %0.3f, %0.3f, %0.3f \n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4),psiParamsFit(5));

end