


run = input('Which run?','s');



%% Initialize Q+

modelType = 'watson';
modelParameters = struct;
modelParameters.tau = 0.5:0.5:2;	% time constant of the center filter (in msecs)
modelParameters.kappa = 0.5:0.5:2;	% multiplier of the time-constant for the surround
modelParameters.zeta = 0:0.5:2;	% multiplier of the amplitude of the surround
modelParameters.beta = 0.5:0.5:2; % multiplier that maps watson 0-1 to BOLD % bins
modelParameters.sigma = 0:0.5:2;	% width of the BOLD fMRI noise against the 0-1 y vals

%{
modelType = 'watson';
modelParameters = struct;
modelParameters.tau = 0.5:0.5:;	% time constant of the center filter (in msecs)
modelParameters.kappa = 0.5:0.25:2;	% multiplier of the time-constant for the surround
modelParameters.zeta = 0:0.25:2;	% multiplier of the amplitude of the surround
modelParameters.beta = 0.5:0.2:2; % multiplier that maps watson 0-1 to BOLD % bins
modelParameters.sigma = 0:0.5:2;	% width of the BOLD fMRI noise against the 0-1 y vals

%}

headroom = .1;
[myQpParams, questDataAtStart] = realTime_Init(modelType,modelParameters,'headroom',headroom);

% Initialize questDataUpdated
questDataUpdated = questDataAtStart;

%{
modelType = 'doe';
modelParameters = struct;
modelParameters.Sr = 0.899:0.025:1.099;
modelParameters.k1 = 0.01:0.005:0.03;
modelParameters.k2 = 0.5:0.05:1;
modelParameters.beta = 0.5:0.1:2; % Amplitude of the scaled response; should converge to unity
modelParameters.sigma = 0:0.1:0.5;	% Standard deviation of the scaled (0-1) noise

headroom = .1;
[myQpParams, questDataAtStart] = realTime_Init(modelType,modelParameters,'headroom',headroom);

%}


% initialize a main data vector
mainDataVectorized = [];

stimDataPath = fullfile(subjectPath,strcat('actualStimuli',run,'.txt'));

% If not at scanner, just load the stimuli we already presented. 
if ~atScanner
    stimData = load(fullfile('/Users','nfuser','Documents','rtQuest','TOME_3021_optimal',horzcat('TOME_3021_run',run),horzcat('stimDataRun',run,'.mat')));
    stimulusVec = stimData.params.stimFreq;
end

mainDataPath = fullfile(subjectPath,'processed',strcat('run',run),strcat('mainDatarun',run));

i = 0;
j = 1;


% Wait for the main data file to exist. 
while ~exist(mainDataPath)
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
    
    % Write new stimulus suggestion for use with play_flash                           
    writeNewStimSuggestion(qpQuery(questDataUpdated),fullfile(subjectPath,'stimLog'));

    j = j + 1;

end

