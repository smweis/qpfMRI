function [psiParamsFit]=doeSimulate(Sr_m, k1_m, k2_m, beta_m, sigma_m, TR, trialLength, qpPres, outNum, seed, varargin)
%% A simulation script to be compiled and run at speed for the DOE temporal model
%
% Syntax:
%  [psiParamsFit]=validate_qpDoETFE_simulate(Sr_m, k1_m, k2_m, beta_m, sigma_m, TR, trialLength, qpPres, outNum)
%
% Description:
%	Takes in simulated parameters for the model, the control (TR, trial
%	length), and whether or not Q+ is in control of things. This is
%	designed to be used with the UF Hipergator code, and so every variable
%	coming in will be interpreted as a string.
%
% Inputs:
%   Sr_m           - doe parameter 1
%   k1_m           - doe parameter 2
%   k2_m           - doe parameter 3
%   beta_m         - Scaling parameter, must be 1.
%   sigma_m        - Noise parameter.
%   TR             - Length of a single TR in milliseconds.
%   trialLength    - Length of a single trial in seconds.
%   qpPres         - 1 for Q+, 0 for random.
%   outNum         - Label for the output CSV file
%   seed           - An integer to use to start the rng seed
% 
% Optional positional arguments:
%   nTrials               - Number of trials
%   stimulusStructDeltaT  - Resolution of the stimulus struct in msecs
%   maxBOLDSimulated      - True size of the max BOLD response
%   maxBOLD               - Initial guess for the max size BOLD response
%   baselineStimulus      - Which stimulus (in freq Hz) is the "baseline"? 
%                           This stimulus should be selected with the 
%                           expectation that the neural response to this
%                           stimulus will be minimal as compared to all 
%                           other stimuli.
%   maxBOLDStimulus       - Which stimulus (in freq Hz) is expected to
%                           return the max value? This is just for the 
%                           initial maxBOLD estimate.
%   nOutcomes             - The number of outcome categories.
%   headroom              - The headroom is the proportion of outcomes that 
%                           are reserved above and below the min and max 
%                           output  to account for noise.
%
% Outputs:
%   psiParamsFit   - 1x5 vector of the psychometric function parameters
%                    returned by the BADS fit.
%
%

% Example: 
%{
[psiParamsFit] = doeSimulate('.98', '.003', '.06', '1.00','.4','800',...,
'12','1','true1','987','30','100','1.6','1.0','0','51','.1');
%}
% 
% 
% TODO: 
% 1. We'd like to vary some of the parameters systematically. It would be
% nice to have a csv that was structured in the format of the input we
% want, then we could call that csv as we are running the simulations.
%
% Something like.... (for example, if we wanted to vary maxBOLDSimulated)
% Header: All variables...
% Row 1:  .98 .003 .06 1.00 .4 800 12 1 false1 30 100 1.6
% Row 2:  .98 .003 .06 1.00 .4 800 12 1 false1 30 100 1.6
% Row 3:  .98 .003 .06 1.00 .4 800 12 1 false1 30 100 1.6
% ...
% Row 30: .98 .003 .06 1.00 .4 800 12 1 false1 30 100 1.6
% Row 31: .98 .003 .06 1.00 .4 800 12 1 false1 30 100 1.5


%% Parse input
p = inputParser;

% Required input
p.addRequired('Sr_m',@ischar);
p.addRequired('k1_m',@ischar);
p.addRequired('k2_m',@ischar);
p.addRequired('beta_m',@ischar);
p.addRequired('sigma_m',@ischar);
p.addRequired('TR',@ischar);
p.addRequired('trialLength',@ischar);
p.addRequired('qpPres',@ischar);
p.addRequired('outNum',@ischar);
p.addRequired('seed',@ischar);

% Replace any defaults then parse the input
% Optional positionalparams
p.addOptional('nTrials','30',@ischar);
p.addOptional('stimulusStructDeltaT','100',@ischar);
p.addOptional('maxBOLDSimulated','1.6',@ischar);
p.addOptional('maxBOLD','1.0',@ischar);
p.addOptional('baselineStimulus','0',@ischar);
p.addOptional('maxBOLDStimulus','30',@ischar);
p.addOptional('nOutcomes','51',@ischar);
p.addOptional('headroom','.1',@ischar);


%% Convert the input parameters 
% Establish qpParams
myQpParams = qpParams;

% Parse
p.parse(Sr_m, k1_m, k2_m, beta_m, sigma_m, TR, trialLength, qpPres,outNum, seed, varargin{:});

nTrials = str2double(p.Results.nTrials); 
stimulusStructDeltaT = str2double(p.Results.stimulusStructDeltaT); 
maxBOLDSimulated = str2double(p.Results.maxBOLDSimulated);
maxBOLD = str2double(p.Results.maxBOLD);
baselineStimulus = str2double(p.Results.baselineStimulus);
maxBOLDStimulus = str2double(p.Results.maxBOLDStimulus);
myQpParams.nOutcomes = str2double(p.Results.nOutcomes);
headroom = str2double(p.Results.headroom);

simulatedPsiParams = [str2double(Sr_m) str2double(k1_m) str2double(k2_m) str2double(beta_m) str2double(sigma_m)]; 
trialLength = str2double(trialLength);
TR = str2double(TR);
seed = str2double(seed);


%% Are we simulating old fashioned constant stimuli or using Q+?

% Check that we have sensible input for the qpPres
assert(str2double(qpPres)==1 || str2double(qpPres)==0,'You used a value for qpPres other than 1 or 0.');

simulateConstantStimuli = logical(str2double(qpPres)); 

%% Some initialization
% Add the stimulus domain. ~Log spaced frequencies between 0 and 30 Hz
myQpParams.stimParamsDomainList = {[baselineStimulus,1.875,3.75,7.5,15,30,60]};

% Create an anonymous function from qpDoETemporalModel in which we
% specify the number of outcomes for the y-axis response
myQpParams.qpPF = @(f,p) qpDoETemporalModel(f,p,myQpParams.nOutcomes,headroom);

% Make the psiParamsDomain
Sr = 0.899:0.025:1.099;
k1 = 0.01:0.04:0.4;
k2 = 0.01:0.04:0.4;
beta = 0.4:0.2:2; % Amplitude of the scaled response; should converge to unity
sigma = 0.3:0.2:1;	% Standard deviation of the scaled (0-1) noise

myQpParams.psiParamsDomainList = {Sr, k1, k2, beta, sigma};

% Beta will converge to 1 as maxBOLD gets closer and closer to the
% simulated maxBOLD. As a result, when simulating data, beta should always
% be set to 1. And, Q+ should always be able to incorporate 1 in its
% domain. Assert these conditions are true. 
assert(simulatedPsiParams(4)==1,'Simulated Beta should always be 1.');
assert(ismember(1,beta),'The domain for beta should always include 1.');


% Pick some random params to simulate if none provided (but set the neural
% noise to zero and beta = 1)
if isempty(simulatedPsiParams)
    simulatedPsiParams = [randsample(Sr,1) randsample(k1,1) randsample(k2,1) 1 0.4];
end

% Derive some lower and upper bounds from the parameter ranges. This is
% used later in maximum likelihood fitting
lowerBounds = [Sr(1) k1(1) k2(1) beta(1) sigma(1)];
upperBounds = [Sr(end) k1(end) k2(end) beta(end) sigma(end)];


% We also want to make sure that the veridical values are actually within
% the domain bounds.
for param = 1:length(simulatedPsiParams)
    assert(lowerBounds(param) < simulatedPsiParams(param) && upperBounds(param) > simulatedPsiParams(param),...
        'Parameter %d is not within the bounds of the parameter domain.',param);
end

% Create and save an rng seed to use for this simulation.
rngSeed = rng(seed);
rngSeed = rng(seed);


% Create a simulated observer with binned output
myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,simulatedPsiParams);

% Initialize Q+
questData = qpInitialize(myQpParams);

% Tack on a continuous output simulated observer to myQpParams
myQpParams.continuousPF = @(f) doeTemporalModel(f,simulatedPsiParams);

% Create a copy of Q+
questDataUntrained = questData;

% Create a stimulusVec to hold the trial across the loops
stimulusVec = nan(1,nTrials);

%% Run simulated trials
for tt = 1:nTrials
    
    % If it is the first three trials we force a baseline or maxBOLD event
    if tt<=3
        if tt == 1 || tt == 3
            stimulusVec(tt) = baselineStimulus;
            fprintf('Initial baseline stimulus: %0.3f\n',stimulusVec(tt));
        else
            stimulusVec(tt) = maxBOLDStimulus;
            fprintf('Initial maxBOLD stimulus: %0.3f\n',stimulusVec(tt));
        end
    else
        if ~simulateConstantStimuli
            % get random stimulus
            stimulusVec(tt) = questData.stimParamsDomain(randi(questData.nStimParamsDomain));
            fprintf('Stimuli chosen randomly: %0.3f\n',stimulusVec(tt));
        else
            % get next stimulus from Q+
            stimulusVec(tt) = qpQuery(questData);
            fprintf('Stimuli chosen by Q+: %0.3f\n',stimulusVec(tt));
        end
    end
    
    % Update maxBOLD with our best guess at the maximum BOLD fMRI response
    % that could be evoked by a stimulus (relative to the baseline
    % stimulus), which is the beta value of the model

    psiParamsIndex = qpListMaxArg(questData.posterior);
    psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
    
    % Only update maxBOLD after we've had at least one maxBOLD trial
    if tt > 2
        maxBOLD = maxBOLD.*psiParamsQuest(4);
        fprintf('Using Q+ fit to generate maxBOLD \nmaxBOLD = %0.3f\n Q+ parameters: %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
        maxBOLD, psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));
    end

    % Create a packet
    thePacket = createPacket('nTrials',tt,...,
        'trialLengthSecs',trialLength,...,
        'stimulusStructDeltaT',stimulusStructDeltaT);

    % Obtain outcomes from tfeUpdate 
    [outcomes] = ...
        tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
        'maxBOLDSimulated',maxBOLDSimulated,...
        'rngSeed',rngSeed.Seed,...,
        'maxBOLD',maxBOLD,...,
        'TRmsecs', TR, ...,
        'noiseSD', simulatedPsiParams(5));

    % Grab a naive copy of questData
    questData = questDataUntrained;

    % Update quest data structure. This is the slow step in the simulation.
    for yy = 1:tt
        questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
    end
       
    
end

%% Adjust Beta/MaxBOLD tradeoff
% For the final parameter estimate, we want to assume that Beta is 1 and
% maxBOLD is whatever it WOULD BE if beta were 1. 

% Grab our current beta estimate is: 
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
betaGuess = psiParamsQuest(4);

% Divide maxBOLD by our beta estimate: (beta / beta) = 1, so
% new maxBOLD = maxBOLD/beta 
maxBOLD = maxBOLD/betaGuess;

% Now run through the fitting steps again with the new maxBOLD
% Create a packet
thePacket = createPacket('nTrials',tt,...,
    'trialLengthSecs',trialLength,...,
    'stimulusStructDeltaT',stimulusStructDeltaT);

% Obtain outcomes from tfeUpdate
[outcomes] = ...
    tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, ...
    'maxBOLDSimulated',maxBOLDSimulated,...
    'rngSeed',rngSeed.Seed,...,
    'maxBOLD',maxBOLD,...,
    'TRmsecs', TR, ...,
    'noiseSD', simulatedPsiParams(5));

% Grab a naive copy of questData
questData = questDataUntrained;

% Update quest data structure. This is the slow step in the simulation.
for yy = 1:tt
    questData = qpUpdate(questData,stimulusVec(yy),outcomes(yy));
end



%% Print some final output to the log

% Find out QUEST+'s estimate of the stimulus parameters, obtained
% on the gridded parameter domain.
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
fprintf('Simulated parameters:              %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
    simulatedPsiParams(1),simulatedPsiParams(2),simulatedPsiParams(3),simulatedPsiParams(4),simulatedPsiParams(5));
fprintf('FINAL Max posterior QUEST+ parameters:   %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4),psiParamsQuest(5));

% Find maximum likelihood fit. Use psiParams from QUEST+ as the starting
% parameter for the search, and impose as parameter bounds the range
% provided to QUEST+.
psiParamsFit = qpFitBads(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
    'lowerBounds', lowerBounds,'upperBounds',upperBounds,...
    'plausibleLowerBounds',lowerBounds,'plausibleUpperBounds',upperBounds);
fprintf('FINAL Maximum likelihood fit parameters: %0.4f, %0.4f, %0.4f, %0.4f, %0.4f \n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4),psiParamsFit(5));
fprintf('FINAL maxBOLD estimate: %0.3f',maxBOLD);

%% Output
T = array2table(psiParamsFit,'VariableNames',{'Sr','k1','k2','beta','sigma'});
T.maxBOLD = maxBOLD;
outfilename = horzcat('doe_',outNum,'.csv');
%save(outfilename,psiParamsFit);
writetable(T,outfilename);

end
