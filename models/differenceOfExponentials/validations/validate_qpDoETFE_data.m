%% QP + DoE TTF + TFE with example data




%% Are we debugging?
reinitializeQuest = 0;

% Clean up
if reinitializeQuest
    clearvars('-except','reinitializeQuest');
    close all;
else
    clearvars('-except','reinitializeQuest','questDataCopy');
    close all;
end



%% Load in fMRI + stim frequency data

% Where is stimulus data?
stimDataLoc = ['processed_fmri_data', filesep, 'optimalResults.mat'];

% Load the data
load(stimDataLoc,'stimParams','detrendTimeseries');

% How many runs?
nRuns = size(detrendTimeseries,1);

% Which stimulus (in freq Hz) is the "baseline" stimulus? This stimulus
% should be selected with the expectation that the neural response to this
% stimulus will be minimal as compared to all other stimuli.
baselineStimulus = 0;




%% Define qpParams

% Some information about the trials?
trialLengthSecs = 12; % seconds per trial
stimulusStructDeltaT = 100; % the resolution of the stimulus struct in msecs
TRmsecs = 800;

% Define the simulated size of the BOLD response and the size that will be
% assumed at the start of the modeling.
fitMaxBOLD = 1.0;

% How talkative is the analysis?
showPlots = true;
verbose = true;


%% Set up Q+

% Get the default Q+ params
myQpParams = qpParams;

% Add the stimulus domain. ~Log spaced frequencies between 0 and 30 Hz
myQpParams.stimParamsDomainList = {unique(stimParams(1).params.stimFreq)};
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
Sr = 0.91:0.05:1.11;
k1 = 0.5:0.05:1;
k2 = 0.005:0.005:0.2;
beta = 0.8:0.05:1.2; % multiplier that maps 0-1 to BOLD % bins
sigma = 1:0.325:2;	% width of the BOLD fMRI noise against the 0-1 y vals
myQpParams.psiParamsDomainList = {Sr, k1, k2, beta, sigma};




%% Create a plot of the entire, integrated dataset
stimulusVecFull = [];
yValsFull = [];
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
    
    [~, ~, ~, yVals] = ...
        tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
        'fitMaxBOLD',fitMaxBOLD);
    
    stimulusVecFull = [stimulusVecFull stimulusVec];
    yValsFull = [yValsFull yVals'];
end

% Obtain the fmincon best fit of the DoE model to the data
myObj = @(p) sqrt(sum((yValsFull-doeTemporalModel(stimulusVecFull,p)).^2));
x0 = [0.8 0.15 0.16 1];
lb = [0 0 0 0];
ub = [3 1 1 3];
unconstrainedFitParams = fmincon(myObj,x0,[],[],[],[],lb,ub);

% Display the data
unconstrainedFit = figure();
stimulusVecPlot = stimulusVecFull;
stimulusVecPlot(stimulusVecPlot==0)=0.01;
freqDomain = logspace(log10(0.01),log10(100),100);
semilogx(freqDomain,doeTemporalModel(freqDomain,unconstrainedFitParams),'-r');
hold on
scatter(stimulusVecPlot,yValsFull,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);

% And the median data values
freqList = unique(stimulusVecPlot);
for ii=1:length(freqList)
    medianVals(ii) = median(yValsFull(stimulusVecPlot==freqList(ii)));
end
plot(freqList,medianVals,'*k')




%% Initialize QP
% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = [Sr(1) k1(1) k2(1) beta(1) sigma(1)];
upperBounds = [Sr(end) k1(end) k2(end) beta(end) sigma(end)];

% Warn the user that we are initializing
if verbose
    tic
    fprintf('Initializing Q+. This may take a minute...\n');
end

% Initialize Q+. Save some time if we're debugging
if ~reinitializeQuest
    if exist('questDataCopy','var')
        questData = questDataCopy;
    else
        questData = qpInitialize(myQpParams);
        questDataCopy = questData;
    end
end

% Prompt the user to start
if verbose
    toc
    fprintf('Press space to start.\n');
    pause
end




%% Loop over runs
maxBOLDPerRun = nan(1,nRuns);
for rr = 1:nRuns

    if verbose
        fprintf('Processing run %d of %d.\n',rr,nRuns);
    end
    
    % Create a copy of Q+
    questDataAtRunStart = questData;
    
    % Obtain the stimulus vec for this run (with zero trials = baseline)
    stimulusVecFull = stimParams(rr).params.stimFreq;
    stimulusVecFull(stimulusVecFull == 0) = baselineStimulus;
    nTrials = length(stimulusVecFull); % how many trials
        
    % Create a plot in which we can track the model progress
    if showPlots
        figure
        
        % Create a full length packet to format the plots
        thePacket = createPacket('nTrials',nTrials,...,
            'trialLengthSecs',trialLengthSecs,...,
            'stimulusStructDeltaT',stimulusStructDeltaT);
        
        % Set up the BOLD fMRI response and model fit
        subplot(3,1,1)
        currentBOLDHandleData = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-k');
        hold on
        currentBOLDHandleFit = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-r');
        xlim([min(thePacket.stimulus.timebase) max(thePacket.stimulus.timebase)]);
        ylim([ floor(min(min(detrendTimeseries))*2)/2, ceil(max(max(detrendTimeseries))*2)/2 ]);
        xlabel('time [msecs]');
        ylabel('BOLD fMRI % change');
        title('BOLD fMRI data');
        
        % Set up the TTF figure
        subplot(3,1,2)
        freqDomain = logspace(log10(0.01),log10(100),100);
        ylim([-0.5 1.5]);
        xlabel('log stimulus Frequency [Hz]');
        ylabel('Relative response amplitude');
        title('Estimate of DoE TTF');
        currentTTFHandle = semilogx(freqDomain,ones(size(freqDomain)),'-k');
        hold on
        currentOutcomesHandle = scatter(nan,nan);
        
        % Calculate the lower headroom bin offset. We'll use this later
        nLower = round(headroom*myQpParams.nOutcomes);
        nUpper = round(headroom*myQpParams.nOutcomes);
        nMid = myQpParams.nOutcomes - nLower - nUpper;
        
        % Set up the entropy x trial figure
        subplot(3,1,3)
        entropyAfterTrial = nan(1,nTrials*nRuns);
        currentEntropyHandle = plot(1:nTrials*nRuns,entropyAfterTrial,'*k');
        xlim([1 nTrials*nRuns]);
        title('Model entropy by trial number');
        xlabel('Trial number');
        ylabel('Entropy');
    end
    
    
    %% Run through trials for this run
    for tt = 1:nTrials
        
        tic
        
        % Obtain the current portion of the full stimulus vec
        stimulusVec = stimulusVecFull(1:tt);
        
        % Skip initial trials until we have a baseline trial available
        foo=1;
        
        % Here is where we would ask QP to supply our next stimulus.
        % qpQuery(questData);
        
        % Update fitMaxBOLD with our best guess at the maximum BOLD
        % fMRI response that could be evoked by a stimulus (relative to the
        % baseline stimulus). The beta value of the model is the 4th parameter.
        % Our hope is that it converges to unity when we have the correct
        % fitMaxBOLD value
        psiParamsIndex = qpListMaxArg(questData.posterior);
        psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
        fitMaxBOLD = fitMaxBOLD.*psiParamsQuest(4);
        
        % Grab a naive copy of questData
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
            'fitMaxBOLD',fitMaxBOLD);
        
        % Update quest data structure. This is the slow step in the simulation.
        for yy = 1:tt
            questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
        end
        
        trialCompTime(rr,tt)=toc;
        
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
            yVals = yVals + mean(yVals(stimulusVec==baselineStimulus));
            stimulusVecPlot = stimulusVec;
            stimulusVecPlot(stimulusVecPlot==0)=0.01;
            delete(currentOutcomesHandle);
            currentOutcomesHandle = scatter(stimulusVecPlot(1:tt),yVals,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
            psiParamsIndex = qpListMaxArg(questData.posterior);
            psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
            delete(currentTTFHandle)
            currentTTFHandle = semilogx(freqDomain,doeTemporalModel(freqDomain,psiParamsQuest),'-r');
            ylim([-0.5 1.5]);

            % Entropy by trial
            subplot(3,1,3)
            delete(currentEntropyHandle)
            entropyAfterTrial=questData.entropyAfterTrial;
            plot(1:length(entropyAfterTrial),entropyAfterTrial,'*k');
            xlim([1 nTrials*nRuns]);
            ylim([0 nanmax(entropyAfterTrial)]);
            xlabel('Trial number');
            ylabel('Entropy');
            
            drawnow
        end % Update plots
        
    end % Loop over trials in this run
    
    maxBOLDPerRun(rr) = fitMaxBOLD;
    
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
fprintf('Max posterior QUEST+ parameters:   %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5),fitMaxBOLD);

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
semilogx(freqDomain,doeTemporalModel(freqDomain,[psiParamsQuest(1:3) mean(maxBOLDPerRun)]),'-k');
semilogx(freqDomain,doeTemporalModel(freqDomain,[psiParamsFit(1:3) mean(maxBOLDPerRun)]),'-m');
legend({'unconstrained fit','data','median vals','Q+ max posterior','Max likelihood'})
