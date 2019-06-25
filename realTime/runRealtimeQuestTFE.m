function [nextStim] = runRealtimeQuestTFE(subject,run,atScanner,model,varargin)

%% Todo! Figure out where the nextStim should be stored. 
%%

% The main function for quest+ TFE data processing in real time at the
% scanner
%
% Syntax:
%   nextStim = runRealtimeQuestTFE(subject,run,atScanner,varargin)
%
% Description:
%	Takes in the subject and run IDs. Churns through the Q+ and TFE
%	pipeline and returns the next stimulus suggestion. 
%
% Inputs:
%   subject               - String. The name/ID of the subject. 
%   run                   - String. The run or acquisition number. Will
%                           generate a folder with the string 'run' before
%                           it. 
%   atScanner             - Logical. Are you actually at the scanner?
%   model                 - String. Name of the model to use to model the
%                           fMRI data ('watson' or 'doe')

% Optional key/value pairs:
%   headroom              - Scalar. Proportion of bins used for headroom by
%                           the Quest+ model. 
% Outputs:
%   nextStim              - Int. Suggestion for the next stimulus parameter
%                           from Quest+



% Examples:

%{

subject = 'Ozzy_Test';
run = '1';
atScanner = true;
headroom = .1;
model = 'watson';
nextStim = runRealtimeQuestTFE(subject,run,atScanner,model,varargin)
% 
%}

%% Parse input
p = inputParser;

% Required input
p.addRequired('subject',@isstr);
p.addRequired('run',@isstr);
p.addRequired('atScanner',@islogical);
p.addRequired('model',@isstr);

% Optional params
p.addParameter('headroom',.1,@isnumeric);

% Parse
p.parse( subject, run, atScanner, model, varargin{:});




%% Get Relevant Paths

[subjectPath, scannerPath, ~, ~] = getPaths(subject);

%% Initialize Q+

if strcmp(model,'watson')

    modelParameters = struct;
    modelParameters.tau = 0.5:0.5:2;	% time constant of the center filter (in msecs)
    modelParameters.kappa = 0.5:0.5:2;	% multiplier of the time-constant for the surround
    modelParameters.zeta = 0:0.5:2;	% multiplier of the amplitude of the surround
    modelParameters.beta = 0.5:0.5:2; % multiplier that maps watson 0-1 to BOLD % bins
    modelParameters.sigma = 0:0.5:2;	% width of the BOLD fMRI noise against the 0-1 y vals

elseif strcmp(model,'doe')
    modelParameters = struct;
    modelParameters.Sr = 0.899:0.025:1.099;
    modelParameters.k1 = 0.01:0.005:0.03;
    modelParameters.k2 = 0.5:0.05:1;
    modelParameters.beta = 0.5:0.1:2; % Amplitude of the scaled response; should converge to unity
    modelParameters.sigma = 0:0.1:0.5;	% Standard deviation of the scaled (0-1) noise
end

headroom = .1;

[myQpParams, questDataAtStart] = realTime_Init(model,modelParameters,'headroom',headroom);



% Initialize questDataUpdated
questDataUpdated = questDataAtStart;


% initialize a main data vector
mainDataVectorized = [];




% CHANGE ME TO THE SCANNER PATH OR WHEREVER ACTUAL STIMULI ARE BEING SAVED! 
stimDataPath = fullfile(subjectPath,strcat('actualStimuli',run,'.txt'));









% If not at scanner, just load the stimuli we already presented. 
if ~atScanner
    stimData = load(fullfile(subjectPath,subject,horzcat('run',run),horzcat('stimDataRun',run,'.mat')));
    stimulusVec = stimData.params.stimFreq;
end




mainDataPath = fullfile(subjectPath,'processed',strcat('run',run),strcat('mainData',run));

i = 0;
j = 1;


% Wait for the main data file to exist. 
while ~exist(mainDataPath,'file')
    continue
    pause(.01);
end
    

while i < 10000000000
    i = i + 1;
    

    load(mainDataPath,'mainData');
   
    
    % Turn main data into a vector
    mainDataVectorized(end+1:end+length(mainData(j).roiSignal)) = mainData(j).roiSignal;

    questDataUpdated = realTime_Update(mainDataVectorized,...
                   stimDataPath,...
                   questDataAtStart, questDataUpdated, myQpParams, modelType,...
                    'headroom',headroom,...
                   'simulationStimVec',stimulusVec);
    
    % CHANGE ME TO THE PATH ON THE SCANNER
    % Write new stimulus suggestion for use with play_flash                           
    writeNewStimSuggestion(qpQuery(questDataUpdated),fullfile(subjectPath,'stimLog'));

    j = j + 1;

end

