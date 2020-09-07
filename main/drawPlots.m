function [mainFig,handleStruct] = drawPlots(myQpfmriParams,myQpParams,stimulusVec,mainFig,handleStruct,thePacketOut,modelResponseStruct,yVals,yValsPlusBaseline,psiParamsQuest,maxBOLDLatestGuess,entropyAfterTrial,varargin)
%% [mainFig,handleStruct] = drawPlots(myQpfmriParams,myQpParams,stimulusVec,mainFig,handleStruct,thePacketOut,modelResponseStruct,yVals,yValsPlusBaseline,psiParamsQuest,maxBOLDLatestGuess,entropyAfterTrial,varargin)

% Draw plots for the main qpfMRI functions
%
% Syntax:
%[mainFig,handleStruct] = drawPlots(myQpfmriParams,myQpParams,stimulusVec,mainFig,handleStruct,thePacketOut,modelResponseStruct,yVals,yValsPlusBaseline,psiParamsQuest,maxBOLDLatestGuess,entropyAfterTrial,varargin)

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
%   mainFig                 - Figure handle. Where to plot. 
%   handleStruct            - Struct. Holder for all plotted features.  
%   thePacketOut            - Struct. Output from tfeUpdate containing the
%                             stimulus struct and simulated response
%                             struct.
%   modelResponseStruct     - Struct. Output from tfeUpdate containing the
%                             fitted response. 
%   yVals                   - Vector. Y-values not normalized to the
%                             baseline stimulus but scaled to maxBOLDLatestGuess. 
%   yValsPlusBaseline       - Vector. Y-values scaled to maxBOLDLatestGuess. 
%   psiParamsQuest          - Vector. Psychometric parameters returned by Q+. 
%   maxBOLDLatestGuess      - Scalar. Latest guess at what maxBOLD is. 
%   entropyAfterTrial       - Vector. Entropy value after each trial. 
% Optional key/value pairs (used in plotting):
%   
% Outputs:
%   mainFig                 - Figure handle. Empty figure to fill and use with drawPlots. 
%   handleStruct            - Struct. All the handles used for plotting for
%                             use with drawPlots.

%% Handle inputs
p = inputParser;

p.addRequired('myQpfmriParams',@isstruct);
p.addRequired('myQpParams',@isstruct);
p.addRequired('stimulusVec',@isvector);
p.addRequired('mainFig',@(x) isa(x,'matlab.ui.Figure'));
p.addRequired('handleStruct',@isstruct);
p.addRequired('yVals',@isvector);
p.addRequired('yValsPlusBaseline',@isvector);
p.addRequired('psiParamsQuest',@isvector);
p.addRequired('maxBOLDLatestGuess',@isscalar);
p.addRequired('entropyAfterTrial',@isvector);
% Optional
p.addParameter('trialTypes',{},@iscell);
p.addParameter('colorStruct',{},@isstruct);

% Optional for saving
p.addParameter('saveGif',false,@islogical);


% Parse inputs
p.parse( myQpfmriParams, myQpParams, stimulusVec, mainFig, handleStruct, yVals, yValsPlusBaseline, psiParamsQuest, maxBOLDLatestGuess, entropyAfterTrial, varargin{:});


%% Setup parameters
% Get nTrials that have been presented
nTrials = length(stimulusVec(~isnan(stimulusVec)));

% Setup colors for the main simulation plot to vary by trial type:
if isempty(p.Results.colorStruct)
    colorStruct = struct;
    colorStruct.baseline = '#0B6F17';
    colorStruct.maxBOLD = '#ADF7B6';
    colorStruct.qPlus = '#007EA7';
    colorStruct.random = '#007EA7';
end

% A cell array to hold cell values as trials come in.
lineColor = cell(1,nTrials);
if isempty(p.Results.trialTypes)
    lineColor{:} = '#007EA7';
else
    for i = 1:nTrials
        if contains(p.Results.trialTypes{i},'baseline')
            lineColor{i} = colorStruct.baseline;
        elseif contains(p.Results.trialTypes{i},'maxbold','IgnoreCase',true)
            lineColor{i} = colorStruct.maxBOLD;
        elseif contains(p.Results.trialTypes{i},'qplus')
            lineColor{i} = colorStruct.qPlus;
        elseif contains(p.Results.trialTypes{i},'random')
            lineColor{i} = colorStruct.random;
        end
    end
end

% Delete all previous annotations and reset the handleStruct
delete(findall(gcf,'type','annotation'));

% Define plot function
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

% Create results strings for anotations
resultString = sprintf('Output value %.03f\n',yVals(end));
if ~isempty(p.Results.trialTypes)
    trialString = sprintf('\nTrial %d: %s Stimulus Selection\n',length(yVals),p.Results.trialTypes{end});
else
    trialString = 'No Trial Type Saved.';
    warning('No Trial Type Data Saved.');
end
stimString = sprintf('Stimulus: %0.3f\n',stimulusVec(length(yVals)));

%% Subplot (left top panel): 
% Veridical model, individual trials, current model fit.
subplot(3,4,[1 2 5 6])

% Plot to the minimum of the stimulus domain excluding zero.
stimulusVecPlot = stimulusVec;
stimulusVecPlot(stimulusVecPlot==0)=min(myQpParams.stimParamsDomainList{1});

% Calculate the predicted relative response from Q+, which
% calculates the model fit depending on the baseline estimate and
% scales to the maxBOLDLatestGuess
predictedQuestRelativeResponse = (myQpfmriParams.model(stimulusDomainFine,psiParamsQuest) - ...
    myQpfmriParams.model(myQpfmriParams.baselineStimulus,psiParamsQuest))*psiParamsQuest(myQpfmriParams.betaIndex)*maxBOLDLatestGuess;

% The current best fit model
delete(handleStruct.currentTTFHandle);
handleStruct.currentTTFHandle = plotFunc(stimulusDomainFine,predictedQuestRelativeResponse,'-r');

% Individual trial outcomes
delete(handleStruct.currentOutcomesHandle);
handleStruct.currentOutcomesHandle = scatter(stimulusVecPlot(1:nTrials-1),...,
    yVals(1:end-1),'o','MarkerFaceColor','b','MarkerEdgeColor','none',...,
    'MarkerFaceAlpha',.2);
% Last individual trial outcome
delete(handleStruct.lastOutcomeHandle);
handleStruct.lastOutcomeHandle = scatter(stimulusVecPlot(nTrials),...,
    yVals(end),'o','MarkerFaceColor','b','MarkerEdgeColor','none',...,
    'MarkerFaceAlpha',1,'HandleVisibility','off');

% Labeling
legend('Veridical Model','Individual Trials','Best-fit Model','Location','northwest');
set(gca,'box','off');

drawnow

%% Subplot (right top panel): 
% Annotations.
subplot(3,4,[3 4]);

% Calculate where to put things
ax = gca;
xLoc = ax.Position(1);
yLoc = ax.Position(2);

% Add annotations
annotation('textbox',[xLoc yLoc .4 .2],'String',sprintf('%s%s%s',trialString,stimString,resultString),'EdgeColor','none','FontSize',15);

% Labeling
title('Trial Information','FontSize',25);

%set(gca,'visible','off');
drawnow

%% Subplot (right middle panel):
% Entropy by trial.
subplot(3,4,[7 8]);

% Plot entropy by trials
handleStruct.currentEntropyHandle=plot(1:nTrials,entropyAfterTrial,'*k');

% Set and keep x and y-lims
xlim([1 myQpfmriParams.nTrials]);
ylim([0 nanmax(entropyAfterTrial)]);

% Labeling
title('Model entropy by trial number');
xlabel('Trial number');
ylabel('Entropy');
set(gca,'box','off');

drawnow

%% Subplot (bottom panel): 
% Simulated fMRI data and model response.
subplot(3,4,[9 10 11 12]);

% In black dots, the simulated TRs. 
delete(handleStruct.currentBOLDHandleData);     
handleStruct.currentBOLDHandleData = plot(thePacketOut.response.timebase./1000,thePacketOut.response.values,'.k');

% In a red line, the model fit. 
delete(handleStruct.currentBOLDHandleFit);
handleStruct.currentBOLDHandleFit = plot(modelResponseStruct.timebase./1000,modelResponseStruct.values,'-r');

% In horizontal bars (colored by trial type), the beta estimate, scaled by
% maxBOLDLatestGuess, per trial.
delete(handleStruct.linePlot);
for m = 1:nTrials
    xLine = [(m-1)*myQpfmriParams.trialLength m*myQpfmriParams.trialLength];
    yLine = [yValsPlusBaseline(m) yValsPlusBaseline(m)];
    handleStruct.linePlot(m) = plot(xLine,yLine,'Color',lineColor{m},'LineWidth',4);
end

% Custom legend handling is the easiest way to go
h = zeros(5,1);
h(1) = plot(NaN,NaN,'.k');
h(2) = plot(NaN,NaN,'-r');
h(3) = plot(NaN,NaN,'-','Color',colorStruct.baseline);
h(4) = plot(NaN,NaN,'-','Color',colorStruct.maxBOLD);
h(5) = plot(NaN,NaN,'-','Color',colorStruct.qPlus);
lgd = legend(h,'TR','Model Fit','Baseline Trial','MaxBOLD Trial','Q+ or Random Selection');
lgd.FontSize = 10;
lgd.NumColumns=2;
legend('boxoff');

drawnow


% Save this plot as a GIF?
if p.Results.saveGif
    if nTrials == 1
        % Requires gif package from fileexchange.
        try
            gif('qpSimulate.gif','DelayTime',1,'frame',gcf,'nodither');
        catch
            warning('You tried to save a gif, but perhaps you are missing the gif package.')
        end
    else
        gif;
    end
end
    