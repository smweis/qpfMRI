function [questDataUpdated, thePacketOut, maxBOLD] = realTime_trialUpdate(timeseries, stimDataPath, questDataAtStart, questDataUpdated, myQpParams, modelType, varargin)
% Takes in the updated timeseries, updates the fits, returns the next trial
% suggestion. Also real-time plots. 
%
% Syntax:
%  [questDataUpdated, thePacketOut] = realTime_Update(timeseries, stimDataPath, questDataAtStart, questDataUpdated, myQpParams, modelType, varargin)
%
% Description:
%	Takes in Q+ objects and parameters and a timeseries of BOLD FMRI data
%   and returns the next trial as suggested by Q+. Also plots the BOLD 
%   timeseries and fit, as well as plots the timeseries, fits, and current best model fit. 
%
% Inputs:
%   timeseries            - 1xn vector of scalars. The BOLD timeseries.
%   stimDataPath          - String. Location where stimulus data will be
%                           saved. 
%   questDataAtStart      - Struct. Initialized questData structure. 
%   questDataUpdated      - Struct. questData structure that has been
%                           updated with all current trials. 
%   myQpParams              - Struct. Generated from myQpParams. This should 
%                           contain a value for nOutcomes other than the
%                           default (2) to ensure enough range of values
%                           for Q+ to work with.
%   modelType             - String. Specify the temporal model used. 
%
% Optional key/value pairs (used in fitting):
%  'baselinestimulus'     - Scalar. The proportion of the nOutcomes from 
%                           myQpParams that will be used as extra on top and
%                           bottom.
%  'maxBOLD'              - Scalar. The value (in % change units) of the
%                           maximum expected response to a stimulus w.r.t.
%                           the response to the baseline stimulus.
%  'headroom'             - Scalar. The proportion of the nOutcomes from 
%                           myQpParams that will be used as extra on top and
%                           bottom.
%
% Optional key/value pairs (used for thePacket, TFE):
%  'TRmsecs'              - Integer. The length of each TR in milliseconds.
%  'trialLengthSecs'      - Integer. The length of each trial in seconds.
%  'stimulusStructDeltaT' - Integer. The resolution of the stimulus struct.
% Optional key/value pairs (for pipeline management)
%   'simulationStimVec'   - If empty, we are in real scanning mode and
%                           actualStimuliRun[run].txt should be used. 
%                           Otherwise, a stimulus vec. 
% 
% Outputs:
%   questDataOut          - Struct. Updated QuestData
%   thePacketOut          - Struct. The updated packet with response struct
%                           completed.
%   nextStim              - Integer. Next stimulus suggested by quest.
% Examples:
%{

%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('timeseries',@isvector);
p.addRequired('stimDataPath', @isstr);
p.addRequired('questDataAtStart',@isstruct);
p.addRequired('questDataUpdated',@isstruct);
p.addRequired('myQpParams',@isstruct);
p.addRequired('modelType',@isstr);

% Optional params used in fitting
p.addParameter('baselinestimulus', 0, @isnumeric);
p.addParameter('maxBOLD', 0.6, @isscalar);
p.addParameter('headroom', .1, @isscalar);

% Optional params used for thePacket (TFE)
p.addParameter('TRmsecs', 800, @isnumeric);
p.addParameter('trialLengthSecs', 12, @isnumeric);
p.addParameter('stimulusStructDeltaT', 100, @isnumeric);


% Optional params used for pipeline management
p.addOptional('simulationStimVec',[]);

% Parse
p.parse(timeseries, stimDataPath, questDataAtStart, questDataUpdated, myQpParams, modelType, varargin{:});



%% Load in the presented stimuli
if isempty(p.Results.simulationStimVec)
    stimulusVec = readActualStimuli(p.Results.stimDataPath);
else
    stimulusVec = p.Results.simulationStimVec;
end

maxBOLD = p.Results.maxBOLD;
baselinestimulus = p.Results.baselinestimulus;
headroom = p.Results.headroom;
TRmsecs = p.Results.TRmsecs;

% Calculate the number of TRs per trial
TRperTrial = p.Results.trialLengthSecs/(TRmsecs/1000);

% How many trials to estimate, based on length of timeseries.
nTrials = ceil(length(timeseries)/TRperTrial);


stimulusVec = stimulusVec(1:nTrials);

% How talkative is the analysis?
showPlots = true;


% Calculate the lower headroom bin offset. We'll use this later
nLower = round(headroom*myQpParams.nOutcomes);
nUpper = round(headroom*myQpParams.nOutcomes);
nMid = myQpParams.nOutcomes - nLower - nUpper;


% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = [myQpParams.psiParamsDomainList{1}(1),...
               myQpParams.psiParamsDomainList{2}(1),...
               myQpParams.psiParamsDomainList{3}(1),...
               myQpParams.psiParamsDomainList{4}(1),...
               myQpParams.psiParamsDomainList{5}(1)];
upperBounds = [myQpParams.psiParamsDomainList{1}(end),...
               myQpParams.psiParamsDomainList{2}(end),...
               myQpParams.psiParamsDomainList{3}(end),...
               myQpParams.psiParamsDomainList{4}(end),...
               myQpParams.psiParamsDomainList{5}(end)];

           

% Create a full length packet to format the plots
thePacket = createPacket('nTrials',nTrials,...,
   'trialLengthSecs',p.Results.trialLengthSecs,...,
   'stimulusStructDeltaT',p.Results.stimulusStructDeltaT);
           
% Create a plot in which we can track the model progress
if showPlots && nTrials > 2
   figure

   % Set up the BOLD fMRI response and model fit
   subplot(3,1,1)
   currentBOLDHandleData = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-k');
   hold on
   currentBOLDHandleFit = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-r');
   xlim([min(thePacket.stimulus.timebase) max(thePacket.stimulus.timebase)]);
   ylim([ floor(min(min(timeseries))*2)/2, ceil(max(max(timeseries))*2)/2 ]);
   xlabel('time [msecs]');
   ylabel('BOLD fMRI % change');
   title('BOLD fMRI data');

   % Set up the TTF figure
   subplot(3,1,2)
   freqDomain = logspace(log10(0.01),log10(100),100);
   currentTTFHandle = semilogx(freqDomain,ones(size(freqDomain)),'-k');
   hold on
   ylim([-0.5 1.5]);
   xlabel('log stimulus Frequency [Hz]');
   ylabel('Relative response amplitude');
   title('Estimate of TTF');
   currentOutcomesHandle = scatter(nan,nan);

   % Set up the entropy x trial figure
   subplot(3,1,3)
   entropyAfterTrial = nan(1,nTrials);
   currentEntropyHandle = plot(1:nTrials,entropyAfterTrial,'*k');
   xlim([1 nTrials]);
   title('Model entropy by trial number');
   xlabel('Trial number');
   ylabel('Entropy');
end


%% Run through trials for this run


tic

% Update maxBOLD with our best guess at the maximum BOLD fMRI
% response that could be evoked by a stimulus (relative to the
% baseline stimulus).
psiParamsIndex = qpListMaxArg(questDataUpdated.posterior);
psiParamsQuest = questDataUpdated.psiParamsDomain(psiParamsIndex,:);
maxBOLD = maxBOLD.*psiParamsQuest(4);



% Add the mean-centered response vector
thePacket.response.values = timeseries;
thePacket.response.values = thePacket.response.values - mean(thePacket.response.values);
thePacket.response.timebase = 0:TRmsecs:length(thePacket.response.values)*TRmsecs - TRmsecs;

% Obtain the outcomes with tfeUpdate
[outcomes, modelResponseStruct, thePacketOut] = ...
    tfeUpdate(thePacket, myQpParams, stimulusVec, baselinestimulus, ...
    'maxBOLD',maxBOLD);

% Update quest data structure. This is the slow step in the simulation.
for yy = 1:nTrials
    questDataAtStart = qpUpdate(questDataAtStart,stimulusVec(yy),outcomes(yy));
end

questDataUpdated = questDataAtStart;

toc

% Update the plots
if showPlots && nTrials > 2
    
    % Simulated BOLD fMRI time-series and fit
    subplot(3,1,1)
    delete(currentBOLDHandleData)
    delete(currentBOLDHandleFit)
    currentBOLDHandleData = plot(thePacketOut.response.timebase,thePacketOut.response.values,'.k');
    currentBOLDHandleFit = plot(modelResponseStruct.timebase,modelResponseStruct.values,'-r');
    
    % TTF figure
    subplot(3,1,2)
    yVals = (outcomes - nLower - 1)./nMid;
    stimulusVecPlot = stimulusVec;
    stimulusVecPlot(stimulusVecPlot==0)=0.01;
    delete(currentOutcomesHandle);
    currentOutcomesHandle = scatter(stimulusVecPlot,yVals,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
    psiParamsIndex = qpListMaxArg(questDataUpdated.posterior);
    psiParamsQuest = questDataUpdated.psiParamsDomain(psiParamsIndex,:);
    delete(currentTTFHandle)
    if strcmpi(modelType,'doe')
        currentTTFHandle = semilogx(freqDomain,doeTemporalModel(freqDomain,psiParamsQuest),'-r');
    elseif strcmpi(modelType,'watson')
        currentTTFHandle = semilogx(freqDomain,watsonTemporalModel(freqDomain,psiParamsQuest),'-r');
    end
    
    ylim([-0.5 1.5]);
    
    % Entropy by trial
    subplot(3,1,3)
    delete(currentEntropyHandle)
    entropyAfterTrial=questDataUpdated.entropyAfterTrial;
    plot(1:length(entropyAfterTrial),entropyAfterTrial,'*k');
    xlim([1 nTrials]);
    ylim([0 nanmax(entropyAfterTrial)]);
    xlabel('Trial number');
    ylabel('Entropy');
    
    drawnow
end % Update plots



%% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questDataUpdated.posterior);
psiParamsQuest = questDataUpdated.psiParamsDomain(psiParamsIndex,:);
fprintf('Max posterior QUEST+ parameters:   %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));

%% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFit(questDataUpdated.trialData,questDataUpdated.qpPF,psiParamsQuest,questDataUpdated.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4),psiParamsFit(5));
