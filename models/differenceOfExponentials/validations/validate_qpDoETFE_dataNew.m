%% QP + DoE TTF + TFE with example data

%% Load in fMRI + stim frequency data

% Where is stimulus data?
stimDataLoc = ['processed_fmri_data', filesep, 'optimalResults.mat'];

% Load the data
load(stimDataLoc,'stimParams','detrendTimeseries');

% Which run?
nRuns = size(detrendTimeseries,1);

% Initial guess for the max size of the evoked BOLD response
maxBOLD = 0.6;

% Which stimulus (in freq Hz) is the "baseline" stimulus? This stimulus
% should be selected with the expectation that the neural response to this
% stimulus will be minimal as compared to all other stimuli.
baselineStimulus = 0;


%% Define experiment specific 

% Some information about the trials?
trialLengthSecs = 12; % seconds per trial
stimulusStructDeltaT = 100; % the resolution of the stimulus struct in msecs
TRmsecs = 800;

% infer nTrials
nTrials = length(detrendTimeseries(1,:))/(trialLengthSecs/TRmsecs*1000);


% How talkative is the analysis?
showPlots = true;
verbose = true;



% Create a packet
thePacket = createPacket('nTrials',nTrials,...,
    'trialLengthSecs',trialLengthSecs,...,
    'stimulusStructDeltaT',stimulusStructDeltaT);


%% Initialize QP

modelType = 'doe';
modelParameters = struct;
modelParameters.Sr = 0.899:0.025:1.099;
modelParameters.k1 = 0.01:0.005:0.03;
modelParameters.k2 = 0.5:0.05:1;
modelParameters.beta = 0.5:0.1:2; % Amplitude of the scaled response; should converge to unity
modelParameters.sigma = 0:0.1:0.5;	% Standard deviation of the scaled (0-1) noise
stimulusDomain = [0,1.875,3.75,7.5,15,20,30];
headroom = .1;

[myQpParams, questData, plottingStruct] = realTime_acquisitionInit(modelType,...
    modelParameters,stimulusDomain,'headroom',headroom,'showPlots',showPlots,'thePacket',thePacket,'nTrials',nTrials);



% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = [modelParameters.Sr(1) modelParameters.k1(1) modelParameters.k2(1) modelParameters.beta(1) modelParameters.sigma(1)];
upperBounds = [modelParameters.Sr(end) modelParameters.k1(end) modelParameters.k2(end) modelParameters.beta(end) modelParameters.sigma(end)];

% Prompt the user to start
if verbose
    fprintf('Press space to start.\n');
    pause
end

% Calculate the lower headroom bin offset. We'll use this later
nLower = round(headroom*myQpParams.nOutcomes);
nUpper = round(headroom*myQpParams.nOutcomes);
nMid = myQpParams.nOutcomes - nLower - nUpper;


%% Create a plot of the entire, integrated dataset
stimulusVecFull = [];
adjustedAmplitudesFull = [];
for rr = 1:nRuns
    stimulusVec = stimParams(rr).params.stimFreq;
    stimulusVec(stimulusVec == 0) = baselineStimulus;
    nTrials = length(stimulusVec); % how many trials
    
    % Create a packet
    thePacket = createPacket('nTrials',nTrials,...,
        'trialLengthSecs',trialLengthSecs,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);
    
    % Add the mean-centered response vector
    thePacket.response.values = detrendTimeseries(rr,1:nTrials*trialLengthSecs/(TRmsecs/1000));
    thePacket.response.values = thePacket.response.values - mean(thePacket.response.values);
    thePacket.response.timebase = 0:TRmsecs:length(thePacket.response.values)*TRmsecs - TRmsecs;
    
    [~, ~, ~, adjustedAmplitudes] = ...
        tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus);

    % Save the values
    adjustedAmplitudesFull = [adjustedAmplitudesFull adjustedAmplitudes'];
    stimulusVecFull = [stimulusVecFull stimulusVec];
end

% Obtain the fmincon best fit of the DoE model to the data
myObj = @(p) sqrt(sum((adjustedAmplitudesFull-doeTemporalModel(stimulusVecFull,p)).^2));
x0 = [0.9 0.1 0.1 1];
unconstrainedFitParams = fminsearch(myObj,x0);

% Display the data
unconstrainedFit = figure;
figure(unconstrainedFit);
stimulusVecPlot = stimulusVecFull;
stimulusVecPlot(stimulusVecPlot==0)=0.01;
freqDomain = logspace(log10(0.01),log10(100),100);
semilogx(freqDomain,doeTemporalModel(freqDomain,unconstrainedFitParams),'-r');
hold on
scatter(stimulusVecPlot,adjustedAmplitudesFull,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);

% And the median data values
freqList = unique(stimulusVecPlot);
for ii=1:length(freqList)
    medianVals(ii) = median(adjustedAmplitudesFull(stimulusVecPlot==freqList(ii)));
end
plot(freqList,medianVals,'*k')
drawnow

hold off;



%% Loop over runs
maxBOLDPerRun = nan(1,nRuns);
for rr = 1:nRuns

    if verbose
        fprintf('Processing run %d of %d.\n',rr,nRuns);
    end
    
    % Save the state of questData at the start of this run
    questDataAtRunStart = questData;
    
    % Obtain the stimulus vec for this run (with zero trials = baseline)
    stimulusVecFull = stimParams(rr).params.stimFreq;
    stimulusVecFull(stimulusVecFull == 0) = baselineStimulus;
    nTrials = length(stimulusVecFull); % how many trials

    
    
    %% Run through trials for this run
    for tt = 1:nTrials
        
        tic
        
        % Obtain the current portion of the full stimulus vec
        stimulusVec = stimulusVecFull(1:tt);
        
        % Skip initial trials until we have a baseline trial available
        % Could add a check here to do so.
        
        % Here is where we would ask QP to supply our next stimulus.
        % qpQuery(questData);
        
        % Update maxBOLD with our best guess at the maximum BOLD fMRI
        % response that could be evoked by a stimulus (relative to the
        % baseline stimulus).
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        maxBOLD = maxBOLD.*psiParamsQuest(4);
        
        % Reset questData to its state at the start of the run
        questData = questDataAtRunStart;
        
        % Create a packet
        thePacket = createPacket('nTrials',tt,...,
            'trialLengthSecs',trialLengthSecs,...,
            'stimulusStructDeltaT',stimulusStructDeltaT);
        
        % Add the mean-centered response vector
        thePacket.response.values = detrendTimeseries(rr,1:tt*trialLengthSecs/(TRmsecs/1000));
        thePacket.response.values = thePacket.response.values - mean(thePacket.response.values);
        thePacket.response.timebase = 0:TRmsecs:length(thePacket.response.values)*TRmsecs - TRmsecs;
        
        % Obtain the outcomes with tfeUpdate
        [outcomes, modelResponseStruct, thePacketOut] = ...
            tfeUpdate(thePacket, myQpParams, stimulusVec(1:tt), baselineStimulus, ...
            'maxBOLD',maxBOLD);
        
        % Update quest data structure. This is the slow step in the simulation.
        for yy = 1:tt
            questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
        end
        
        trialCompTime(rr,tt)=toc;
        
        % Update the plots
        if showPlots
            
            figure(plottingStruct.figureHandle);
            
            % Simulated BOLD fMRI time-series and fit
            subplot(3,1,1)
            delete(plottingStruct.currentBOLDHandleData);
            delete(plottingStruct.currentBOLDHandleFit);
            plottingStruct.currentBOLDHandleData = plot(thePacketOut.response.timebase,thePacketOut.response.values,'.k');
            plottingStruct.currentBOLDHandleFit = plot(modelResponseStruct.timebase,modelResponseStruct.values,'-r');
            
            % TTF figure
            subplot(3,1,2)
            yVals = (outcomes - nLower - 1)./nMid;
            stimulusVecPlot = stimulusVec;
            stimulusVecPlot(stimulusVecPlot==0)=0.01;
            delete(plottingStruct.currentOutcomesHandle);
            plottingStruct.currentOutcomesHandle = scatter(stimulusVecPlot(1:tt),yVals,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
            psiParamsIndex = qpListMaxArg(questData.posterior);
            psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
            delete(plottingStruct.currentTTFHandle);
            plottingStruct.currentTTFHandle = semilogx(plottingStruct.freqDomain,doeTemporalModel(plottingStruct.freqDomain,psiParamsQuest),'-r');
            ylim([-0.5 1.5]);

            % Entropy by trial
            subplot(3,1,3)
            delete(plottingStruct.currentEntropyHandle);
            plottingStruct.entropyAfterTrial=questData.entropyAfterTrial;
            plottingStruct.currentEntropyHandle = plot(1:length(plottingStruct.entropyAfterTrial),plottingStruct.entropyAfterTrial,'*k');
            xlim([1 nTrials*nRuns]);
            ylim([0 nanmax(plottingStruct.entropyAfterTrial)]);
            xlabel('Trial number');
            ylabel('Entropy');
            
            drawnow
        end % Update plots
        
    end % Loop over trials in this run
    
    maxBOLDPerRun(rr) = maxBOLD;
    
end % Loop over runs

% Done with the simulation
if verbose
    fprintf('\n');
    fprintf('Median trial update time: %2.2f seconds\n',median(median(trialCompTime)));
end

% Plot the computation time per trial
if showPlots
    figure
    plot(median(trialCompTime),'xr');
    xlabel('trial number');
    ylabel('computation time [secs]');
end


%% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
fprintf('Max posterior QUEST+ parameters:   %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));

%% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFit(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4),psiParamsFit(5));

%% Update the unconstrained fit plot
% Display the data
figure(unconstrainedFit);
freqDomain = logspace(log10(0.01),log10(100),100);
semilogx(freqDomain,doeTemporalModel(freqDomain,[psiParamsQuest(1:3) psiParamsQuest(4)*mean(maxBOLDPerRun)]),'-k');
semilogx(freqDomain,doeTemporalModel(freqDomain,[psiParamsFit(1:3) psiParamsFit(4)*mean(maxBOLDPerRun)]),'-m');
legend({'unconstrained fit','data','median vals','Q+ max posterior','Max likelihood'})
