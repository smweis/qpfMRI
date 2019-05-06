function [myQpParams, questDataAtStart, plottingStruct] = realTime_acquisitionInit(modelType,modelParameters,stimulusDomain,varargin)
% Initialize Q+ for temporal fitting engine and real-time use. Can also set
% up figures. 
%
% Syntax:
%  [myQpParams, questDataAtStart, plottingStruct] = realTime_acquisitionInit(modelType,modelParameters,stimulusDomain,varargin)

%
% Description:
%	Initializes a Q+ object for use with Watson or DOE modeling of realtime
%	fMRI data.
%
% Inputs:
%  modelType              - String. The name of the model used in fitting.
%                           For now either doe or watson.
%  modelParameters        - Struct. Wherein each field (named with that 
%                           parameter name) is a vector of the possible 
%                           values Q+ could return for that parameter.
%  stimulusDomain         - Vector. 1xk vector of integers where k is the number of
%                           possible stimulus values that could be presented. 


%
% Optional key/value pairs (Q+ parameters):
%  'nOutcomes'            - Integer. Number of possible bins Q+ could use.

%  'headroom'             - Scalar. The proportion of the nOutcomes from 
%                           qpParams that will be used as extra on top and
%                           bottom.
%  'showPlots'            - Logical. If True, pass out handles to
%                           initalized figures.
%  'thePacket'            - TFE packet. Necessary if showing plots.
%  'nTrials'              - Integer. Number of trials to be presented.
%                           Necessary if showing plots. 
%
% Outputs:
%   myQpParams            - Struct. Parameter settings for use to initalize
%                           Q+. We may need some of these later.
%   questDataAtStart      - Initialized struct for use with Q+
%   plottingStruct        - Struct. 
%                           3x1 Figure handle and subplots. 
%                           subplot 1: timeseriesFigure - Plots the fMRI 
%                               timeseries data and model fit from TFE.
%                           subplot 2: modelFigure - Plots the stimulus 
%                               domain by BOLD response along with the 
%                               temporal function fit. 
%                           subplot 3: entropyFigure - Plots the time course
%                               of the entropy in temporal function model space. 
%
% Examples:
%{
% Using the DOE model
modelType = 'doe';
modelParameters = struct;
modelParameters.Sr = 0.899:0.025:1.099;
modelParameters.k1 = 0.01:0.005:0.03;
modelParameters.k2 = 0.5:0.05:1;
modelParameters.beta = 0.5:0.1:2; % Amplitude of the scaled response; should converge to unity
modelParameters.sigma = 0:0.1:0.5;	% Standard deviation of the scaled (0-1) noise

stimulusDomain = [0,1.875,3.75,7.5,15,20,30];


[myQpParams, questDataAtStart, plottingStruct] = realTime_acquisitionInit(modelType,modelParameters,stimulusDomain);


% Using the Watson model
modelType = 'watson';
modelParameters = struct;
modelParameters.tau = 0.5:0.5:8;	% time constant of the center filter (in msecs)
modelParameters.kappa = 0.5:0.25:2;	% multiplier of the time-constant for the surround
modelParameters.zeta = 0:0.25:2;	% multiplier of the amplitude of the surround
modelParameters.beta = 0.5:0.2:2; % multiplier that maps watson 0-1 to BOLD % bins
modelParameters.sigma = 0:0.5:2;	% width of the BOLD fMRI noise against the 0-1 y vals

stimulusDomain = [0,1.875,3.75,7.5,15,20,30];

[myQpParams, questDataAtStart, plottingStruct] = realTime_acquisitionInit(modelType,modelParameters,stimulusDomain);

%}


%% Parse input
p = inputParser;

p.addRequired('modelType',@isstr);
p.addRequired('modelParameters',@isstruct);
p.addRequired('stimulusDomain',@isvector);

% Optional params used for Q+
p.addParameter('nOutcomes', 51, @isnumeric);
p.addParameter('headroom', 0.1, @isnumeric);
p.addParameter('showPlots', 1, @islogical);
p.addParameter('thePacket',[]);
p.addParameter('nTrials', 30 ,@isnumeric);


% Parse
p.parse(modelType,modelParameters,stimulusDomain,varargin{:});





%% Set up Q+

% Get the default Q+ params
myQpParams = qpParams;

% Add the stimulus domain.
myQpParams.stimParamsDomainList = {stimulusDomain};

% The number of outcome categories.
myQpParams.nOutcomes = p.Results.nOutcomes;


headroom = p.Results.headroom;

% Create an anonymous function depending on model type in which we
% specify the number of outcomes for the y-axis response
if strcmpi(modelType,'doe')
    myQpParams.qpPF = @(f,p) qpDoETemporalModel(f,p,myQpParams.nOutcomes,headroom);
    paramNames = {'Sr','k1','k2','beta','sigma'};
elseif strcmpi(modelType,'watson')
    myQpParams.qpPF = @(f,p) qpWatsonTemporalModel(f,p,myQpParams.nOutcomes,headroom);
    paramNames = {'tau','kappa','zeta','beta','sigma'};
else
    error('The model specification you have entered has not been created yet. Please try \n DOE or Watson.');
end

myQpParams.psiParamsDomainList = {modelParameters.(paramNames{1}),...
                                  modelParameters.(paramNames{2}),...
                                  modelParameters.(paramNames{3}),...
                                  modelParameters.(paramNames{4}),...
                                  modelParameters.(paramNames{5})};


%% Initialize QP and create a packet


tic
fprintf('Initializing Q+. This may take a minute...\n');


questDataAtStart = qpInitialize(myQpParams);
        

toc



% Set up the BOLD fMRI response and model fit
if p.Results.showPlots
    
    thePacket = p.Results.thePacket;
    
    plottingStruct = struct;

    plottingStruct.figureHandle = figure; 

    % Set up the BOLD fMRI response and model fit
    subplot(3,1,1)
    plottingStruct.currentBOLDHandleData = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-k');
    hold on
    plottingStruct.currentBOLDHandleFit = plot(thePacket.stimulus.timebase,zeros(size(thePacket.stimulus.timebase)),'-r');
    xlim([min(thePacket.stimulus.timebase) max(thePacket.stimulus.timebase)]);
    ylim([-2 2]);
    xlabel('time [msecs]');
    ylabel('BOLD fMRI % change');
    title('BOLD fMRI data');

    % Set up the TTF figure
    subplot(3,1,2)
    plottingStruct.freqDomain = logspace(log10(0.01),log10(100),100);
    plottingStruct.currentTTFHandle = semilogx(plottingStruct.freqDomain,ones(size(plottingStruct.freqDomain)),'-k');
    hold on
    ylim([-0.5 1.5]);
    xlabel('log stimulus Frequency [Hz]');
    ylabel('Relative response amplitude');
    title('Estimate of DoE TTF');
    plottingStruct.currentOutcomesHandle = scatter(nan,nan);

    % Set up the entropy x trial figure
    subplot(3,1,3)
    plottingStruct.entropyAfterTrial = nan(1,p.Results.nTrials);
    plottingStruct.currentEntropyHandle = plot(1:p.Results.nTrials,plottingStruct.entropyAfterTrial,'*k');
    hold on
    xlim([1 p.Results.nTrials]);
    title('Model entropy by trial number');
    xlabel('Trial number');
    ylabel('Entropy');

end
    

