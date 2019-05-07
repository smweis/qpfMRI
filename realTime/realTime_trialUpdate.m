function [questDataOut, thePacketOut, nextStim, maxBOLDOut, plottingStruct] = realTime_trialUpdate(timeseries, stimDataPath, questDataAtStart, myQpParams, modelFunction, plottingStruct, varargin)
% Takes in the updated timeseries, updates the fits, returns the next trial
% suggestion. Also real-time plots. 
%
% Syntax:
%  [questDataUpdated, thePacketOut, nextStim, maxBOLD, plottingStruct] = realTime_trialUpdate(timeseries, stimDataPath, questDataAtStart, myQpParams, modelFunction, plottingStruct, varargin)
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
%   myQpParams              - Struct. Generated from myQpParams. This should 
%                           contain a value for nOutcomes other than the
%                           default (2) to ensure enough range of values
%                           for Q+ to work with.
%   modelFunction          - Function handle. Specify the temporal model used. 
%   plottingStruct         - Struct. 
%                           3x1 Figure handle and subplots. 
%                           subplot 1: timeseriesFigure - Plots the fMRI 
%                               timeseries data and model fit from TFE.
%                           subplot 2: modelFigure - Plots the stimulus 
%                               domain by BOLD response along with the 
%                               temporal function fit. 
%                           subplot 3: entropyFigure - Plots the time course
%                               of the entropy in temporal function model space. 
%
% Optional key/value pairs (used in fitting):
%  'baselineStimulus'     - Scalar. The proportion of the nOutcomes from 
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
%
% Optional key/value pairs (for pipeline management)
%   'simulationStimVec'   - If empty, we are in real scanning mode and
%                           actualStimuliRun[run].txt should be used. 
%                           Otherwise, a stimulus vec. 
%   'showPlots'            - Logical. If True, pass out handles to
%                           initalized figures.

% 
% Outputs:
%   questDataOut          - Struct. Updated QuestData
%   thePacketOut          - Struct. The updated packet with response struct
%                           completed.
%   nextStim              - Integer. Next stimulus suggested by quest.
%   maxBOLD               - Scalar. Estimate for maximum BOLD response.
%   plottingStruct        - Struct. 
%                           3x1 Figure handle and subplots. 
%                           subplot 1: timeseriesFigure - Plots the fMRI 
%                               timeseries data and model fit from TFE.
%                           subplot 2: modelFigure - Plots the stimulus 
%                               domain by BOLD response along with the 
%                               temporal function fit. 
%                           subplot 3: entropyFigure - Plots the time course
%                               of the entropy in temporal function model space. 
% Examples:
%{

%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('timeseries',@isvector);
p.addRequired('stimulusVec', @isvector);
p.addRequired('questDataAtStart',@isstruct);
p.addRequired('myQpParams',@isstruct);
p.addRequired('modelFunction',@(x) isa(x,'function_handle'));
p.addRequired('plottingStruct',@isstruct);

% Optional params used in fitting
p.addParameter('baselineStimulus', 0, @isnumeric);
p.addParameter('maxBOLD', 0.6, @isscalar);
p.addParameter('headroom', .1, @isscalar);

% Optional params used for thePacket (TFE)
p.addParameter('TRmsecs', 800, @isnumeric);
p.addParameter('trialLengthSecs', 12, @isnumeric);
p.addParameter('stimulusStructDeltaT', 100, @isnumeric);

% Optional params, other
p.addParameter('showPlots', 1, @islogical);

% Parse
p.parse(timeseries, stimDataPath, questDataAtStart, myQpParams, modelFunction, plottingStruct, varargin{:});



%% Load in the presented stimuli


% Calculate the number of TRs per trial
TRperTrial = p.Results.trialLengthSecs/(p.Results.TRmsecs/1000);

% How many trials to estimate, based on length of timeseries. Trim stimvec
% to that length. 
nTrials = ceil(length(timeseries)/TRperTrial);
stimulusVec = p.Results.stimulusVec(1:nTrials);


% Calculate the lower headroom bin offset. We'll use this later
nLower = round(p.Results.headroom*myQpParams.nOutcomes);
nUpper = round(p.Results.headroom*myQpParams.nOutcomes);
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


%% Run through trials for this run

% Initialize questDataOut
questDataOut = questDataAtStart; 

tic

% Add the mean-centered response vector
thePacket.response.values = timeseries;
thePacket.response.values = thePacket.response.values - mean(thePacket.response.values);
thePacket.response.timebase = 0:p.Results.TRmsecs:length(thePacket.response.values)*p.Results.TRmsecs - p.Results.TRmsecs;

% Obtain the outcomes with tfeUpdate
[outcomes, modelResponseStruct, thePacketOut] = ...
    tfeUpdate(thePacket, myQpParams, stimulusVec, p.Results.baselineStimulus, ...
    'maxBOLD',p.Results.maxBOLD);

% Update quest data structure. This is the slow step in the simulation.
for yy = 1:nTrials
    questDataOut = qpUpdate(questDataOut,stimulusVec(yy),outcomes(yy));
end



% Update maxBOLD with our best guess at the maximum BOLD fMRI
% response that could be evoked by a stimulus (relative to the
% baseline stimulus).
psiParamsIndex = qpListMaxArg(questDataOut.posterior);
psiParamsQuest = questDataOut.psiParamsDomain(psiParamsIndex,:);
maxBOLDOut = p.Results.maxBOLD.*psiParamsQuest(4);

nextStim = qpQuery(questDataOut);


toc

% Update the plots
if p.Results.showPlots && nTrials > 2

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
    plottingStruct.currentOutcomesHandle = scatter(stimulusVecPlot(1:nTrials),yVals,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
    psiParamsIndex = qpListMaxArg(questDataOut.posterior);
    psiParamsQuest = questDataOut.psiParamsDomain(psiParamsIndex,:);
    delete(plottingStruct.currentTTFHandle);
    plottingStruct.currentTTFHandle = semilogx(plottingStruct.freqDomain,doeTemporalModel(plottingStruct.freqDomain,psiParamsQuest),'-r');
    ylim([-0.5 1.5]);

    % Entropy by trial
    subplot(3,1,3)
    delete(plottingStruct.currentEntropyHandle);
    plottingStruct.entropyAfterTrial=questDataOut.entropyAfterTrial;
    plottingStruct.currentEntropyHandle = plot(1:length(plottingStruct.entropyAfterTrial),plottingStruct.entropyAfterTrial,'*k');
    xlim([1 nTrials]);
    ylim([0 nanmax(plottingStruct.entropyAfterTrial)]);
    xlabel('Trial number');
    ylabel('Entropy');

    drawnow
end % Update plots



%% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questDataOut.posterior);
psiParamsQuest = questDataOut.psiParamsDomain(psiParamsIndex,:);
fprintf('Max posterior QUEST+ parameters:   %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));

%% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFit(questDataOut.trialData,questDataOut.qpPF,psiParamsQuest,questDataOut.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4),psiParamsFit(5));
