function [mainFig,handleStruct] = initializePlots(myQpfmriParams,myQpParams,varargin)
%% [mainFig,handleStruct] = initializePlots(myQpfmriParams,myQpParams,varargin)
% Initialize plots for the main qpfMRI functions
%
% Syntax:
%  [mainFig,handleStruct] = initializePlots(myQpfmriParams,myQpParams,varargin)
%
% Description:
%	Takes in a model and a possible set of parameters and whether or not Q+ 
%   is in control of things flag. 
%
% Inputs:
% Required Inputs
%   myQpfmriParams          - Struct. Set of parameters used for qpfmri.
%                             See qpfmriParams function for more details
%   myQpParams              - Struct. Set of parameters used for qPlus.
%                             See qpParams function for more details
% Optional key/value pairs (used in plotting):
%   'figWidth'              - Integer (Default = 900)
%                             Width of figure window size.
%   'figHeight'             - Integer (Default = 900)
%                             Height of figure window size.
%   'realData'              - Logical. (Default = false)
%                             If true, will not print simulated parameters.
% Outputs:
%   mainFig                 - Figure handle. Empty figure to fill and use with drawPlots. 
%   handleStruct            - Struct. All the handles used for plotting for
%                             use with drawPlots.

p = inputParser;

p.addRequired('myQpfmriParams',@isstruct);
p.addRequired('myQpParams',@isstruct);
p.addParameter('figWidth',1100,@isnumeric);
p.addParameter('figHeight',1100,@isnumeric);
p.addParameter('realData',false,@islogical);

% Parse inputs
p.parse( myQpfmriParams, myQpParams, varargin{:});
   
   
   
if strcmpi(myQpfmriParams.stimulusDomainSpacing,'log')
    plotFunc = @semilogx;
else
    plotFunc = @plot;
end

% Create a fine version of the stimulus space.
if strcmpi(myQpfmriParams.stimulusDomainSpacing,'log')
    stimulusDomainFine = logspace(log(min(myQpParams.stimParamsDomainList{1})),...,
        log(max(myQpParams.stimParamsDomainList{1})),100);
else
    stimulusDomainFine = linspace(min(myQpParams.stimParamsDomainList{1}),...,
        max(myQpParams.stimParamsDomainList{1}),100);
end

% Create an empty packet for plotting
thePacket = createPacket(myQpfmriParams,myQpfmriParams.nTrials);

% Initialize the main figure
mainFig = figure('Position',[10 10 p.Results.figWidth p.Results.figHeight]);
set(gcf,'color','w');
hold on;

%% Subplot (bottom panel): Simulated fMRI data and model response.
subplot(3,4,[9 10 11 12])
handleStruct.currentBOLDHandleData = plot(thePacket.stimulus.timebase./1000,zeros(size(thePacket.stimulus.timebase)),'-k');
hold on
handleStruct.currentBOLDHandleFit = plot(thePacket.stimulus.timebase./1000,zeros(size(thePacket.stimulus.timebase)),'-r');
xlim([min(thePacket.stimulus.timebase./1000) max(thePacket.stimulus.timebase)./1000]);
% Use yMax to set the range for the plot.
yMax = myQpfmriParams.maxBOLDSimulated + myQpfmriParams.maxBOLDSimulated*myQpfmriParams.simulatedPsiParams(myQpfmriParams.sigmaIndex)*2;
ylim([-yMax yMax]);
xlabel('time [seconds]','FontSize',20);
ylabel('BOLD fMRI [% change]','FontSize',20);
title('Simulated BOLD fMRI data','FontSize',20);
set(gca,'box','off');
ax = gca;
ax.FontSize = 15;
drawnow

%% Subplot (left top panel): Veridical model, individual trials, current model fit.
subplot(3,4,[1 2 5 6])
set(gca,'Box','off');
% TO DO: MIGHT WANT TO FIX THIS DOWN THE ROAD??
% Predicted relative response is adjusted for the baseline and scaled
% to maxBOLD
if ~p.Results.realData
    predictedRelativeResponse = (myQpfmriParams.model(stimulusDomainFine,myQpfmriParams.simulatedPsiParams) - ...
        myQpfmriParams.model(myQpfmriParams.baselineStimulus,myQpfmriParams.simulatedPsiParams))*myQpfmriParams.maxBOLDSimulated;
    plotFunc(stimulusDomainFine,predictedRelativeResponse,'-k');


end

ylim([-0.5 myQpfmriParams.maxBOLDSimulated+1]);

xlabel('Stimulus Values');
ylabel('Relative response amplitude');
title('Estimate of Model');
hold on
% Scatter plot of the stimulus values
handleStruct.currentOutcomesHandle = scatter(nan,nan);
handleStruct.lastOutcomeHandle = scatter(nan,nan);
% Veridical model plot.
handleStruct.currentTTFHandle = plot(stimulusDomainFine,myQpfmriParams.model(stimulusDomainFine,myQpfmriParams.simulatedPsiParams),'-k');
handleStruct.linePlot = plot(nan);


%% Subplot (right top panel): Annotations.
subplot(3,4,[3 4]);
title('Trial Information','FontSize',25);
set(gca,'visible','off');
drawnow

%% Subplot (right middle panel): Entropy by trial.
subplot(3,4,[7 8]);
entropyAfterTrial = nan(1,myQpfmriParams.nTrials);
xlabel('Trial number');
ylabel('Entropy');
handleStruct.currentEntropyHandle = plot(1:myQpfmriParams.nTrials,entropyAfterTrial,'*k');
xlim([1 myQpfmriParams.nTrials]);
set(gca,'box','off');
drawnow
end

