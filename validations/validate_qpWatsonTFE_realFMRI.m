%% QP + Watson TTF + TFE + real fmri data


%% Do we want to reinitialize?
reinitializeQuest = 1;

% Clean up
if reinitializeQuest
    clearvars('-except','questDataCopy','reinitializeQuest');
    close all;
else
    clearvars('-except','reinitializeQuest');
    close all;
end

debugFlag = 1;



%% Load in fMRI + stim frequency data and create tfeObj and thePacket


% Let's do this with 1 run for now. 
runNum = 2; % possible values: ints 1-5

% Where is stimulus data?
% Optimal results
stimDataLoc = ['processed_fmri_data', filesep, 'optimalResults.mat'];

% Real-time pipline results
% stimDataLoc = ['processed_fmri_data', filesep, 'rtSimResults.mat'];


load(stimDataLoc);

%% Define the veridical model params

% Leave the simulatedPsiParams empty to try a random set of params.
% Here are some params to try for specific TTF shapes:
%  A low-pass TTF in noisy fMRI data: [10 1 0.83 1]
%  A band-pass TTF in noisy fMRI data: [1.47 1.75 0.83 1]
%simulatedPsiParams = [4 1 1 1 0];
simulatedPsiParams = [];



% Which stimulus (in freq Hz) is the "baseline" stimulus? This stimulus
% should be selected with the expectation that the neural response to this
% stimulus will be minimal as compared to all other stimuli.
baselineStimulus = 200;



% Create the stimulus vec (with no zero trials)
stimulusVec = stimParams(runNum).params.stimFreq;
stimulusVec(stimulusVec == 0) = baselineStimulus;

% Some information about the trials?
nTrials = length(stimulusVec); % how many trials
trialLengthSecs = 12; % seconds per trial
baselineTrialRate = 6; % present a gray screen (baseline trial) every X trials
stimulusStructDeltaT = 100; % the resolution of the stimulus struct in msecs
boldTRmsecs = 800; % msecs


% Use this for calculating the stim struct resolution
stimulusStructPerTrial = trialLengthSecs * 1000 / stimulusStructDeltaT;
responseStructPerTrial = trialLengthSecs * 1000 / boldTRmsecs;

% Initialize thePacket
thePacket = createPacket('nTrials',nTrials,...,
                         'trialLengthSecs',trialLengthSecs,...,
                         'stimulusStructDeltaT',stimulusStructDeltaT);


% Set up the response struct;
thePacket.response.values = detrendTimeseries(runNum,:);
responseTimebaseLength = length(thePacket.stimulus.timebase)/(boldTRmsecs/stimulusStructDeltaT);
thePacket.response.timebase = linspace(0,(responseTimebaseLength-1)*boldTRmsecs,responseTimebaseLength);

% We'll see how quickly we can converge on the full model
watsonParams = watsonParams(runNum,:);
watsonParams(:,4) = 1;

% How talkative is the simulation?
showPlots = true;
verbose = true;

options = optimoptions(@fmincon,'Display', 'off');


%% Set up Q+

% Get the default Q+ params
myQpParams = qpParams;

% Add the stimulus domain. Log spaced frequencies between ~2 and 30 Hz
myQpParams.stimParamsDomainList = {[stimParams(1).params.allFreqs]};
nStims = length(myQpParams.stimParamsDomainList{1}); 

% The number of outcome categories.
myQpParams.nOutcomes = 51;

% The headroom is the proportion of outcomes that are reserved above and
% below the min and max output of the Watson model to account for noise
headroom = .1;

% Create an anonymous function from qpWatsonTemporalModel in which we
% specify the number of outcomes for the y-axis response
myQpParams.qpPF = @(f,p) qpWatsonTemporalModel(f,p,myQpParams.nOutcomes,headroom);

% Define the parameter ranges
tau = 0.5:0.5:8;	% time constant of the center filter (in msecs)
kappa = 0.5:0.25:2;	% multiplier of the time-constant for the surround
zeta = 0:0.25:2;	% multiplier of the amplitude of the surround
beta = 0.8:0.1:1; % multiplier that maps watson 0-1 to BOLD % bins
sigma = 0:0.5:2;	% width of the BOLD fMRI noise against the 0-1 y vals
myQpParams.psiParamsDomainList = {tau, kappa, zeta, beta, sigma};


% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = [tau(1) kappa(1) zeta(1) beta(1) sigma(1)];
upperBounds = [tau(end) kappa(end) zeta(end) beta(end) sigma(end)];

% Create a simulated observer
myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,watsonParams);

% Warn the user that we are initializing
if verbose
    tic
    fprintf('Initializing Q+. This may take a minute...\n');
end

% Initialize Q+. Save some time if we're debugging

if reinitializeQuest
    if exist('questDataCopy','var')
        questData = questDataCopy;
    else
        questData = qpInitialize(myQpParams);
        questDataCopy = questData;
    end
else
    questData = qpInitialize(myQpParams);
    questDataCopy = questData;
end




% Prompt the user we to start the simulation
if verbose
    toc
    fprintf('Press space to start.\n');
    pause
    fprintf('Fitting...');
end

% Create a plot in which we can track the model progress
if showPlots
    figure

    % Set up the BOLD fMRI response and model fit
    subplot(2,1,1)
    currentBOLDHandleData = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-k');
    hold on
    currentBOLDHandleFit = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-r');
    xlim([min(thePacket.stimulus.timebase) max(thePacket.stimulus.timebase)]);
    ylim([-2 3]);    
    xlabel('time [msecs]');
    ylabel('BOLD fMRI % change');
    title('BOLD fMRI data');
    
    % Set up the TTF figure
    subplot(2,1,2)
    freqDomain = logspace(0,log10(100),100);
    semilogx(freqDomain,watsonTemporalModel(freqDomain,watsonParams),'-k');
    ylim([-0.5 1.5]);
    xlabel('log stimulus Frequency [Hz]');
    ylabel('Relative response amplitude');
    title('Estimate of Watson TTF');
    hold on
    currentFuncHandle = semilogx(freqDomain,watsonTemporalModel(freqDomain,watsonParams),'-k');

    % Calculate the lower headroom bin offset. We'll use this later
    nLower = round(headroom*myQpParams.nOutcomes);
    nUpper = round(headroom*myQpParams.nOutcomes);
    nMid = myQpParams.nOutcomes - nLower - nUpper;
    
    % Calculate how many non baseline trials there will be. 
    nNonBaselineTrials = nTrials - ceil(nTrials/baselineTrialRate);


end



nonBaselineTrials = 0;
thePacketOrig = thePacket;

% Create a copy of Q+
questDataUntrained = questData;

% Initiate BOLD limits fairly large
BOLDminFit = -1;
BOLDmaxFit = 1;

BOLDminSim = min(thePacket.response.values) - std(thePacket.response.values);
BOLDmaxSim = max(thePacket.response.values) + std(thePacket.response.values);


%% Run simulated trials
for tt = 1:nTrials
    
    
    % Get stimulus for this trial
    questStim(tt) = qpQuery(questData);
    
    
    questData = questDataUntrained;
    
    stim(tt) = stimulusVec(tt);
    
    % Update thePacket to be just the current trials.
    thePacket.stimulus.values = thePacketOrig.stimulus.values(1:tt,1:tt*stimulusStructPerTrial);
    thePacket.stimulus.timebase = thePacketOrig.stimulus.timebase(:,1:tt*stimulusStructPerTrial);
    
    thePacket.response.values = thePacketOrig.response.values(1:tt*responseStructPerTrial);
    thePacket.response.timebase = thePacketOrig.response.timebase(:,1:tt*responseStructPerTrial);
    
    % Simulate outcome with tfe
    [outcomes, modelResponseStruct, thePacketOut]  = tfeUpdate(thePacket,myQpParams,stim,200);
    
    
    BOLDminFit = min(thePacket.response.values)
    BOLDmaxFit = max(thePacket.response.values)
    
    
    % Update quest data structure
    
    for yy = 1:nonBaselineTrials
        questData = qpUpdate(questData,stim(yy),outcome(yy));
    end
    
    
    
    yVals = (outcomes - nLower - 1)./nMid;
    yVals = yVals + mean(yVals(stim==baselineStimulus));
    stimulusVecPlot = stim;
    stimulusVecPlot(stimulusVecPlot==0)=1;
    
    % This is a lot of code to try to fit the data to obtain watson
    % parameters for plotting. We end up fitting the average of each
    % stimulus frequency across all presentations of that frequency.
    
    % Only works with more than 1 trial of data
    if tt > 2
        
        % Identify the unique stims and the mean of the BOLD response for those
        % stims
        [uniqueStims,~,k] = unique(stim(1,:));
        numberUniqueStims = numel(uniqueStims);
        meanBoldPerStim = zeros(size(stim,1),numberUniqueStims);
        for nu = 1:numberUniqueStims
            indexToThisUniqueValue = (nu==k)';
            meanBoldPerStim(:,nu) = mean(yVals(indexToThisUniqueValue));
            stdBoldPerStim(:,nu) = std(yVals(indexToThisUniqueValue));
        end
        
        stimulusFreqHzFine = logspace(log10(min(stim)),log10(max(stim)),100);
        splineInterpolatedMax = max(spline(uniqueStims,meanBoldPerStim,stimulusFreqHzFine));
        % Scale the x vector so that the max is zero
        meanBoldPerStim = meanBoldPerStim ./ splineInterpolatedMax;
        yVals = yVals ./splineInterpolatedMax;
        myObj = @(p) sqrt(sum((meanBoldPerStim-watsonTemporalModel(uniqueStims,p)).^2));
        x0 = [2 2 2 1];
        watsonParams = fmincon(myObj,x0,[],[],[],[],[],[],[], options);
        
    
    
    end
    





    
    % Update the plot
    if tt > 2 && showPlots
        
        % Current guess at the TTF, along with stims and outcomes

        % Simulated BOLD fMRI time-series and fit       
        subplot(2,1,1)
        delete(currentBOLDHandleData)
        delete(currentBOLDHandleFit)
        currentBOLDHandleData = plot(thePacket.response.timebase,thePacket.response.values,'.k');
        currentBOLDHandleFit = plot(modelResponseStruct.timebase,modelResponseStruct.values,'-r');
        
        
        % TTF figure
        subplot(2,1,2)
        scatter(stim(tt),yVals(tt),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2)
        delete(currentFuncHandle)
        currentFuncHandle = plot(freqDomain,watsonTemporalModel(freqDomain,watsonParams),'-r');

        drawnow
    end
    
    

end
    
% Done with the simulation
if verbose
    fprintf('\n');
end

%{
%% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
fprintf('Simulated parameters: %0.1f, %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    simulatedPsiParams(1),simulatedPsiParams(2),simulatedPsiParams(3),simulatedPsiParams(4),simulatedPsiParams(5));
fprintf('Max posterior QUEST+ parameters: %0.1f, %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));

%% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFit(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4),psiParamsFit(5));

%}


