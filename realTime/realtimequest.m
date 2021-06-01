function [qpfmriResults]=realtimequest(varargin)
%% [qpfmriResults]=realTime_Quest(myQpfmriParams, varargin)
% A script that will take in a BOLD timeseries and output a stimulus
% suggestion
%
% Syntax:
%  [qpfmriResults]=realtimequest(varargin)
%
% Description:
%	Takes in a model and a possible set of parameters and whether or not Q+ 
%   is in control of things flag. 
%
% Inputs:
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

subject = 'Ozzy_Test';
% Provide a model handle
model = @logistic;

% Specify the parameter domain. Each value must correspond to a parameter
% expected by the model. 

paramsDomain = struct;
paramsDomain.slope = makeDomain(-1.2,-.2,10,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.beta = makeDomain(.75,1.25,11,'spacing','zeno');

% Sigma in the parameter domain is searching for noiseSD
paramsDomain.sigma = makeDomain(0.01,1,10);

% Specify a stimulus domain and whether it spaced linear or log.
stimulusDomain = {makeDomain(.5,1,12)};
stimulusDomainSpacing = 'log';

% Number of trials to run and TR.
nTrials = 20;
TRmsecs = 800;

% Allow Q+ to control the stimuli or not (false).
qpPres = true;

% Set the number of outcome categories / bins.
nOutcomes = 15;

% Do you want to see plots?
showPlots = true; 

%How long the trials are (in seconds).
trialLength = 12;

outNum = 1;

[myQpfmriParams,myQpParams] = qpfmriParams(model,paramsDomain,'qpPres',qpPres,...,
'stimulusDomain',stimulusDomain,'stimulusDomainSpacing',stimulusDomainSpacing,...,
'nTrials',nTrials,'trialLength',trialLength,'nOutcomes',nOutcomes,'TR',TRmsecs,'outNum',outNum);


% Run the simulation. 
[qpfmriResults]=realTime_Quest(myQpfmriParams,myQpParams,subject,'showPlots',showPlots);

%}


%% TO DO 
% 
% 2. I'd like to clean up how RANDOM vs. Q+ is handled. 
%       Ideally for random, it's not doing much...maybe fitting the model
%       and plotting just to keep things consistent. 

%% Parse Inputs
debug = 0;
p = inputParser;

% p.addRequired('myQpfmriParams',@isstruct);
% p.addRequired('myQpParams',@isstruct);
% p.addRequired('subject',@isstr);
% 
% Optional params for plotting
p.addParameter('showPlots',false,@islogical);
p.addParameter('saveFigs',false,@islogical);
% % Parse inputs
% p.parse( myQpfmriParams, myQpParams, subject, varargin{:});


if debug
    p.addParameter('myQpfmriParams',@isstruct);
    p.addParameter('myQpParams',@isstruct);
    % Parse inputs
    p.parse(varargin{:});
    myQpfmriParams = p.Results.myQpfmriParams;
    myQpParams = p.Results.myQpParams;
else
    [subject,myQpfmriParams,myQpParams] = qpgetparams();
    % Parse inputs
    p.parse(varargin{:});
end

run = num2str(myQpfmriParams.outNum);
%% Double check a few key parameters

function validateinput(inputVar)
    % Take a string that matches the name of a variable. Tell the user what
    % that variable is currently entered as. They can validate or change.
    
    fprintf('%s is %d.\nCorrect? Press Enter.\n',inputVar,myQpfmriParams.(inputVar));
    pOutput = input('Wrong? Input the correct number: ','s');
    
    if isempty(pOutput)
        fprintf('%s set to %d\n',inputVar, myQpfmriParams.(inputVar));
    elseif floor(str2double(pOutput)) == ceil(str2double(pOutput))
	    myQpfmriParams.(inputVar) = str2double(pOutput);
        fprintf('%s has CHANGED: New value is %d\n',inputVar,myQpfmriParams.(inputVar));
    else
        error("Not an integer");
    end
    fprintf('\n');
end

validateinput('outNum');
validateinput('nTrials');
validateinput('trialLength');
validateinput('TR');

[~, ~, ~, ~, ~, subjectPath] = getpaths(subject, 'neurofeedback');

% dataPath = input("Enter the full path of the directory location of the TIMESERIES: ",'s');
% stimPath = input("Enter the full path of the directory location where the
% STIMULUS SUGGESTION should be written: ",'s');
dataPath = fullfile(subjectPath,'processed',strcat('run',run),'fMRITimeseries_0.txt');
stimPath = fullfile(subjectPath,'stims',strcat('run',run));
stimFile = 'suggestions.txt';
presentedStimFile = strcat('actualStimuli',run,'.txt');
% TODO: add warning if stimFile exists (will be overwritten)


%% Initialize struct to save results out.
qpfmriResults = myQpfmriParams;

% Initialize maxBOLD
qpfmriResults.maxBOLDoverTrials = nan(1,myQpfmriParams.nTrials);
maxBOLDLatestGuess = myQpfmriParams.maxBOLDInitialGuess;
qpfmriResults.psiParamsQuest = nan(myQpfmriParams.nTrials,length(fieldnames(myQpfmriParams.paramsDomain)));
qpfmriResults.entropyOverTrials = cell(1,myQpfmriParams.nTrials);

% Set up save info and directory
folderName = ['.' filesep myQpfmriParams.outFolder];
fileName = ['fMRIresults_' num2str(myQpfmriParams.outNum) '.mat'];
modelFile = ['model_' num2str(myQpfmriParams.outNum) '.mat'];
if ~exist(folderName,'dir')
    mkdir(folderName);
end


%% Some checking on the parameters and making sure Q+ has the models it needs.

% Derive some lower and upper bounds from the parameter ranges. This is
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
lowerBoundsConstrained(myQpfmriParams.betaIndex) = .99999999;
upperBoundsConstrained(myQpfmriParams.betaIndex) = 1.0000001;

% Make sure 1 is a member of the beta parameter domains. If it's not, add
% it in. 
if ~ismember(1,myQpfmriParams.paramsDomain.beta)
    warning('The domain for beta should always include 1. Adding 1 to beta parameter domain. This may create non-linear spacing.');
    myQpfmriParams.paramsDomain.beta = sort(myQpfmriParams.paramsDomain.beta);
    myQpfmriParams.paramsDomain.beta = [myQpfmriParams.paramsDomain.beta(myQpfmriParams.paramsDomain.beta<1) 1 myQpfmriParams.paramsDomain.beta(myQpfmriParams.paramsDomain.beta>1)];
end

% Veridical values should be within domain bounds, except for sigma.
% for param = 1:length(myQpfmriParams.paramNamesInOrder)
%     if ~strcmp(myQpfmriParams.paramNamesInOrder{param},'sigma')
%         assert(lowerBounds(param) <= myQpfmriParams.simulatedPsiParams(param)...
%             && upperBounds(param) >= myQpfmriParams.simulatedPsiParams(param),...,
%             'Parameter %s is not within the bounds of the parameter domain.',myQpfmriParams.paramNamesInOrder{param});
%     end
% end
% 
% Create a simulated observer with binned output
% myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,myQpfmriParams.simulatedPsiParams);

% Calculate the lower headroom bin offset. We'll use this for plotting.
nLower = max([1 round(myQpfmriParams.headroom*myQpParams.nOutcomes)]);
nUpper = max([1 round(myQpfmriParams.headroom*myQpParams.nOutcomes)]);
nMid = myQpParams.nOutcomes - nLower - nUpper;

%% Initialize Q+

% Do we already have data?
% These try-catch blocks will search for .mat files named with the previous
% run number. 
loadToggle = input("Do you want to load previous Q+ data? [Y/n]: ",'s');
if strcmp(loadToggle,'y')
    try
        % Load in the last run's questData
        oldModelFile = ['model_' num2str(myQpfmriParams.outNum - 1) '.mat'];
        load(fullfile(folderName,oldModelFile),'questData');
        fprintf('Re-loading data from run %d\n',myQpfmriParams.outNum);

        % Clean out the the trial-specific stuff. Keep the rest. 
        questData.entropyAfterTrial = [];
        questData.trialData = [];
        questData.stimIndices = [];

    catch
        % Initialize Q+. 
        % Warn the user that initializing has started and may take a minute.
        fprintf('No model file found: %s\n',oldModelFile);
        tic
        fprintf('Initializing Q+. This may take a minute...\n');
        questData = qpInitialize(myQpParams);
        toc
    end

    try
        % Load in last run's maxBOLDFinal estimate.
        oldmaxBold = ['fMRIresults_' num2str(myQpfmriParams.outNum - 1) '.mat'];
        maxBOLDLastRun = load(fullfile(folderName,oldmaxBold),'qpfmriResults');
        maxBOLDLatestGuess = maxBOLDLastRun.qpfmriResults.maxBOLDFinal;
    catch
        fprintf('No results file found: %s\n',oldmaxBold);
        fprintf('Starting with maxBOLD = %d.\n',qpfmriResults.maxBOLDInitialGuess);
    end
else
    tic
    fprintf('Initializing Q+. This may take a minute...\n');
    questData = qpInitialize(myQpParams);
    toc
    fprintf('Starting with maxBOLD = %d.\n',qpfmriResults.maxBOLDInitialGuess);
end

% Tack on a continuous output simulated observer to myQpParams
% myQpParams.continuousPF = @(f) myQpfmriParams.model(f,myQpfmriParams.simulatedPsiParams);

% Store this untrained copy. We will re-train the run model each time.
questDataUntrained = questData;

%% Wait for the initial trials to be completed. 
% Once they're done, store the values for use in the main loop below.
disp("Press ENTER to begin");
pause;

disp(horzcat('Waiting for first ',num2str(myQpfmriParams.baselineMaxBOLDInitial),' trials'));
initTrials = myQpfmriParams.baselineMaxBOLDInitial;

% Create a stimulusVec to hold the trial across the loops
qpfmriResults.stimulusVec = nan(1,initTrials);
% This allows us to store how each trial was generated
qpfmriResults.stimulusVecTrialTypes = cell(1,initTrials);

%% Plotting features
% Create a plot in which we can track the model progress
if p.Results.showPlots
    [mainFig,handleStruct] = initializePlots(myQpfmriParams,myQpParams,'realData',true);
end

%% MAIN LOOP
for iTrial = 1:myQpfmriParams.nTrials
    
    % read in presented stimuli
%     try
%         qpfmriResults.stimulusVec = readmatrix(fullfile(stimPath,presentedStimFile));
%     catch
%         warning('Unable to read in presented stimuli');
%     end
    
    % First trials require baseline or maxBOLD. 
    if iTrial <= myQpfmriParams.baselineMaxBOLDInitial
        % Baseline trial
        if mod(iTrial,2) > 0 
            qpfmriResults.stimulusVec(iTrial) = myQpfmriParams.baselineStimulus;
            qpfmriResults.stimulusVecTrialTypes{iTrial} = 'baseline';
        % maxBOLD trial
        else 
            qpfmriResults.stimulusVec(iTrial) = myQpfmriParams.maxBOLDStimulus;
            qpfmriResults.stimulusVecTrialTypes{iTrial} = 'maxBOLD';
        end
    % Optionally, every X trials can alternate as baseline or maxBOLD. 
    % Baseline trial
    elseif mod(iTrial,myQpfmriParams.baselineMaxBOLDRepeating*2) == 0
        qpfmriResults.stimulusVec(iTrial) = myQpfmriParams.baselineStimulus;
        qpfmriResults.stimulusVecTrialTypes{iTrial} = 'baseline';
    % maxBOLD trial
    elseif mod(iTrial,myQpfmriParams.baselineMaxBOLDRepeating) == 0
        qpfmriResults.stimulusVec(iTrial) = myQpfmriParams.maxBOLDStimulus;
        qpfmriResults.stimulusVecTrialTypes{iTrial} = 'maxBOLD';
    % Every other trial will be selected randomly, or by Q+
    else
        if ~myQpfmriParams.qpPres
            qpfmriResults.stimulusVec(iTrial) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            qpfmriResults.stimulusVecTrialTypes{iTrial} = 'random';
        else
            qpfmriResults.stimulusVec(iTrial) = qpQuery(questData);
            qpfmriResults.stimulusVecTrialTypes{iTrial} = 'qplus';
        end
    end
    
    disp(['Trial #',num2str(iTrial),': ',num2str(qpfmriResults.stimulusVec(iTrial)),...,
        ' ', qpfmriResults.stimulusVecTrialTypes{iTrial}]);
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if iTrial > 2 && ~isnan(qpfmriResults.psiParamsQuest(iTrial-1,myQpfmriParams.betaIndex))
        maxBOLDLatestGuess = maxBOLDLatestGuess.*qpfmriResults.psiParamsQuest(iTrial-1,myQpfmriParams.betaIndex);
    end

    % Create a packet
    thePacket = createPacket(myQpfmriParams,iTrial);
    
    % Load in the timeseries
    if iTrial > 2
        % NOTE: both the timeseries and stimulusVec are trimmed to exclude
        % partial trial data (tt-1 excludes the current trial)
        if myQpfmriParams.qpPres
            try
                timeseries = readmatrix(dataPath);
                % Trim it based on the # of TRs seen in full trials.
                timeseries = timeseries(1:(iTrial-1)*(myQpfmriParams.trialLength*1000/myQpfmriParams.TR));
                %timeseries = timeseries(1:iTrial-1);
                
                
                thePacket.response.values = timeseries;
                thePacket.response.timebase = 0:myQpfmriParams.TR:(length(thePacket.response.values)-1)*myQpfmriParams.TR;
            catch e
                disp(e.message);
                warning("Could not load timeseries");
                if iTrial > myQpfmriParams.baselineMaxBOLDInitial
                    warning("Defaulting to random stimulus");
                    qpfmriResults.stimulusVec(iTrial) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
                    qpfmriResults.stimulusVecTrialTypes{iTrial} = 'random';
                end
            end
        end
        
        % Obtain outcomes from tfeUpdate
        [outcomes, modelResponseStruct, thePacketOut, ~, baselineEstimate] = ...
            tfeUpdate(thePacket, myQpParams, myQpfmriParams, ...,
            qpfmriResults.stimulusVec(1:iTrial-1), maxBOLDLatestGuess);
        
        % Grab a naive copy of questData
        questData = questDataUntrained;
        
        % Update quest data structure. This is the slow step in the simulation.
        for yy = 1:iTrial-1
            questData = qpUpdate(questData,qpfmriResults.stimulusVec(yy),outcomes(yy));
        end
        
        % Save individual run output.
        psiParamsIndex = qpListMaxArg(questData.posterior);
        qpfmriResults.psiParamsQuest(iTrial,:) = questData.psiParamsDomain(psiParamsIndex,:);
        qpfmriResults.maxBOLDoverTrials(iTrial) = maxBOLDLatestGuess;
        qpfmriResults.entropyOverTrials{iTrial} = questData.entropyAfterTrial;
        
        %% Calculate and print results
        % Calculate BOLD values from outcomes, scaled to maxBOLD and the
        % baseline estimate.
        yVals = (outcomes - nLower - 1)./nMid;
        yVals = yVals*qpfmriResults.psiParamsQuest(iTrial,myQpfmriParams.betaIndex)*maxBOLDLatestGuess;
        yValsPlusBaseline = yVals + baselineEstimate;
        resultString = sprintf('Output value %.03f\n',yVals(end));
        fprintf(resultString);
    end
    
    % Print results to the console.
    trialString = sprintf('\nTrial %d: %s Stimulus Selection\n',iTrial,qpfmriResults.stimulusVecTrialTypes{iTrial});
    fprintf(trialString);
    stimString = sprintf('Stimulus: %0.3f\n',qpfmriResults.stimulusVec(iTrial));
    fprintf(stimString);

    fid = fopen(fullfile(stimPath,stimFile),'at');
    fprintf(fid,'%0.3f\n',qpfmriResults.stimulusVec(iTrial));
    fclose(fid);
    
    %% Plot the ongoing results
    % Update the plots
    if p.Results.showPlots
        % Draw plots
        [mainFig,handleStruct] = drawPlots(myQpfmriParams,myQpParams,...,
            qpfmriResults,mainFig,handleStruct,thePacketOut,modelResponseStruct,...,
            yVals,yValsPlusBaseline);
    end
    pause(myQpfmriParams.trialLength);

    
end

%% Adjust Beta/MaxBOLD tradeoff

%% Do the fit
% Grab our current beta estimate: 
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,qpfmriResults.psiParamsQuest(iTrial,:),questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds,...
    'plausibleLowerBounds',lowerBounds,'plausibleUpperBounds',upperBounds);
maxBOLDLatestGuess = maxBOLDLatestGuess.*psiParamsFit(myQpfmriParams.betaIndex);

% Now run through the fitting steps again with the new maxBOLD
thePacket = createPacket(myQpfmriParams,myQpfmriParams.nTrials);
timeseries = readmatrix(dataPath);
timeseries = timeseries(1:(iTrial-1)*(myQpfmriParams.trialLength*1000/myQpfmriParams.TR));
thePacket.response.values = timeseries;
thePacket.response.timebase = 0:myQpfmriParams.TR:(length(thePacket.response.values)-1)*myQpfmriParams.TR;
[outcomes, ~, ~, ~, ~] = ...
        tfeUpdate(thePacket, myQpParams, myQpfmriParams, ...,
        qpfmriResults.stimulusVec, maxBOLDLatestGuess);
    
        
questData = questDataUntrained;
for yy = 1:iTrial
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

% Print Q+ Fit parameters.
fprintf('\nMax posterior QUEST+ parameters:   '); 
for i = 1:length(qpfmriResults.psiParamsQuestFinal)
    fprintf('%s: %0.3f ',myQpfmriParams.paramNamesInOrder{i},qpfmriResults.psiParamsQuestFinal(i));
end

% Prepare for BADS fit.
% For the final parameter estimate, we want to assume that Beta is 1 and
% maxBOLD is whatever it WOULD BE if beta were 1. 
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

% Save out the stim results and the quest struct.
save(fullfile(folderName,fileName),'qpfmriResults');
save(fullfile(folderName,modelFile),'questData');


end