function [outcomes, modelResponseStruct, thePacket, adjustedAmplitudes, baselineEstimate] = tfeUpdate(thePacket, qpParams, myQpfmriParams, stimulusVec, maxBOLDLatestGuess, varargin)
%% [outcomes, modelResponseStruct, thePacket, adjustedAmplitudes, baselineEstimate] = tfeUpdate(thePacket, qpParams, myQpfmriParams, stimulusVec, maxBOLDLatestGuess, varargin)
%Returns the QP outcomes given a packet and a stimulus vector.
%
% Syntax:
%  [outcomes, modelResponseStruct, thePacket, adjustedAmplitudes] = tfeUpdate(thePacket, qpParams, stimulusVec, baselineStimulus)
%
% Description:
%	Takes in the tfeObject created with tfeInit along with thePacket. If
%	thePacket.response is empty, will simulate an fMRI signal, fit that
%	signal, and return outputs suitable for use with Quest +.
%
% Inputs:
%   thePacket             - Struct. Describes the stimulus, observed
%                           response, and (optionally) a convolution kernel
%                           to be applied to the stimulus. If the response
%                           field is empty, the routine will run in
%                           simulation mode.
%   qpParams              - Struct. Generated from qpParams. This should 
%                           contain a value for nOutcomes other than the
%                           default (2) to ensure enough range of values
%                           for Q+ to work with.
%   myQpfmriParams        - Struct. Set of parameters used for qpfmri.
%                           See qpfmriParams function for more details
%   stimulusVec           - 1xn vector. Stimulus values presented thus far.
%   maxBOLDLatestGuess    - Scalar. The latest guess for the value of
%                           maxBOLD
%
% Optional key/value pairs:
%  'rngSeed'              - Struct. By passing a seed to the random
%                           number generator, calling function can ensure
%                           that this routine returns the same output for a
%                           given input in simulation model.
%  'noiseSD'              - Scalar. The amplitude of the noise added to the
%                           simulated BOLD fMRI signal in units of standard
%                           deviations.
%  'pinkNoise'            - Logical. If the noise should be modeled as
%                           pink (auto-correlated).
%  'definePreScanState'   - Logical. Extends the first trial backwards into
%                           negative time (before the start of the scan).
%                           Without this, the model doesn't know how to
%                           treat the neural state prior to the start of
%                           the scan.
%
% Outputs:
%   outcomes              - 1xk vector. The outcome of each stimulus,
%                           expressed as a integer indicating into which of
%                           the nOutcome "bins" the response to each
%                           stimulus fell.
%   modelResponseStruct   - Struct. The simulated response struct from tfe 
%                           method fitResponse
%   thePacket             - Struct. The updated packet with response struct
%                           completed.
%   adjustedAmplitudes    - 
%   baselineEstimate      - 
% Examples:
%{
% Provide a model handle
model = @logistic;

% Specify the parameter domain. Each value must correspond to a parameter
% expected by the model. 

paramsDomain = struct;
paramsDomain.slope = makeDomain(-1.2,-.2,10,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.beta = makeDomain(.75,1.25,11,'spacing','zeno');

% Sigma in the parameter domain is searching for noiseSD
paramsDomain.sigma = makeDomain(.5,4,8);

[myQpfmriParams,myQpParams] = qpfmriParams(model,paramsDomain);


% SIMULATION MODE
thePacket = createPacket(myQpfmriParams,myQpfmriParams.nTrials);

% Generate a random stimulus vector
stimulusDomain = {makeDomain(.01,1,25)};
stimulusVec = randsample(stimulusDomain{:},myQpfmriParams.nTrials,true);
% Need at least one baseline.
stimulusVec(1) = myQpfmriParams.baselineStimulus;

%% Tack on a continuous output simulated observer to myQpParams
myQpParams.continuousPF = @(f) myQpfmriParams.model(f,myQpfmriParams.simulatedPsiParams);

% Perform the simulation
[outcomes, modelResponseStruct, thePacketOut, ~, ~] = tfeUpdate(thePacket,...,
 myQpParams, myQpfmriParams, stimulusVec,...,
 myQpfmriParams.maxBOLDInitialGuess);

% Plot the results
figure
subplot(2,1,1)
plot(thePacketOut.response.timebase,thePacketOut.response.values,'.k');
hold on
plot(modelResponseStruct.timebase,modelResponseStruct.values,'-r');
subplot(2,1,2)
stimulusVecPlot = stimulusVec;
stimulusVecPlot(stimulusVecPlot==0)=1;
semilogx(stimulusVecPlot,outcomes,'xk')
refline(0,myQpParams.nOutcomes.*headroom);
refline(0,myQpParams.nOutcomes-myQpParams.nOutcomes.*myQpfmriParams.headroom);
ylim([1 myQpParams.nOutcomes]);

%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('thePacket',@isstruct);
p.addRequired('qpParams',@isstruct);
p.addRequired('myQpfmriParams',@isstruct);
p.addRequired('stimulusVec',@isvector);
p.addRequired('maxBOLDLatestGuess',@isnumeric);

% Optional params used in simulation
p.addParameter('rngSeed',rng(1,'twister'),@isstruct);
p.addParameter('pinkNoise',0, @isnumeric);
p.addParameter('definePreScanState',true, @islogical);

% Parse
p.parse( thePacket, qpParams, myQpfmriParams, stimulusVec, maxBOLDLatestGuess, varargin{:});

% We need to have at least one "baseline" stimulus in the vector of
% stimuli to support reference-coding of the BOLD fMRI responses
if isempty(find(stimulusVec==myQpfmriParams.baselineStimulus, 1))
    error('The stimulusVec must have at least one instance of the baselineStimulus.');
end

% Setting this within tfeUpdate is to ensure the seed is set properly. 
rngSeed = rng(p.Results.rngSeed);
rngSeed = rng(p.Results.rngSeed);

% Construct the temporal fitting engine model object
tfeObj = tfeIAMP('verbosity','none');

% Set a default params value based on how many stimulus values there should
% have been (which is based on the number of rows in the stimulus.values
% struct)
defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);


%% Handle pre-scan time
% If the stimulus state is not set prior to the start of the scan (i.e., at
% negative time points), extend the state at the earliest time point back
% in time. This extension back is set to be equal to the duration of the
% kernel, so that after convolution the system is in a steady-state prior
% to the scan start.
if min(thePacket.stimulus.timebase)>=0 && p.Results.definePreScanState
    
    %% Stimulus
    
    % Define the pre-scan timebase
    deltaT = diff(thePacket.stimulus.timebase(1:2));
    preScanStart = -ceil(max(thePacket.kernel.timebase)/deltaT)*deltaT;
    preScanEnd = min(thePacket.stimulus.timebase)-deltaT;
    preScan.timebase = preScanStart:deltaT:preScanEnd;
    
    % Extend the stimulus
    preScan.stimValues = repmat(thePacket.stimulus.values(:,1),1,length(preScan.timebase));
    thePacket.stimulus.timebase = [preScan.timebase thePacket.stimulus.timebase];
    thePacket.stimulus.values = [preScan.stimValues thePacket.stimulus.values];
    
    %% Response
    
    %%%%%%%%%%%%%%%%%
    %%
    %% I haven't had the chance to test this with an actual passed response
    %% so if this routine is crashing when presented with empirical data,
    %% chck here for bugs
    %%
    %%%%%%%%%%%%%%%%%
    if ~isempty(thePacket.response)
        
        % Define the pre-scan timebase
        deltaT = diff(thePacket.response.timebase(1:2));
        preScanStart = -ceil(max(thePacket.kernel.timebase)/deltaT)*deltaT;
        preScanEnd = min(thePacket.response.timebase)-deltaT;
        preScan.timebase = preScanStart:deltaT:preScanEnd;
        
        % Extend the response
        preScan.respValues = repmat(thePacket.response.values(:,1),1,length(preScan.timebase));
        thePacket.response.timebase = [preScan.timebase thePacket.response.timebase];
        thePacket.response.values = [preScan.respValues thePacket.response.values];
    end
end


%% If response is empty, simulate it
if isempty(thePacket.response)
    
    % Initialize params0, which will allow us to create the forward model.
    params0 = tfeObj.defaultParams('defaultParamsInfo', defaultParamsInfo);
    params0.noiseSd = myQpfmriParams.fMRInoise;
    params0.noiseInverseFrequencyPower = p.Results.pinkNoise;
    modelAmplitude = zeros(length(stimulusVec),1);
        
    % Obtain the continuous amplitude response
    for ii = 1:length(stimulusVec)        
        modelAmplitude(ii) = qpParams.continuousPF(stimulusVec(ii)) .* myQpfmriParams.maxBOLDSimulated;    
    end

    % We enforce reference coding, such that amplitude of response to the
    % stimuli is expressed relative to the response to the baseline
    % stimulus.
    modelAmplitude = modelAmplitude - mean(modelAmplitude(stimulusVec==myQpfmriParams.baselineStimulus));

    % Place the responses in the paramMainMatrix
    params0.paramMainMatrix = modelAmplitude;
        
    % Lock the MATLAB random number generator to give us the same BOLD
    % noise on every iteration.
    rng(p.Results.rngSeed);
    
    % Create a simulated BOLD fMRI time series
    thePacket.response = tfeObj.computeResponse(params0,thePacket.stimulus,thePacket.kernel,'AddNoise',true);
    
    % Resample the response to the BOLD fMRI TR
    BOLDtimebase = min(thePacket.response.timebase):myQpfmriParams.TR:max(thePacket.response.timebase);
    thePacket.response= tfeObj.resampleTimebase(thePacket.response,BOLDtimebase);
    
    % Mean center the response
    thePacket.response.values = thePacket.response.values - mean(thePacket.response.values);
    
end


%% Fit the response
[params,~,modelResponseStruct] = tfeObj.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');


% We engage in reference coding, such that the amplitude of any stimulus is
% expressed w.r.t. the "baseline" stimulus
baselineEstimate = mean(params.paramMainMatrix(stimulusVec==myQpfmriParams.baselineStimulus));
adjustedAmplitudes = params.paramMainMatrix - baselineEstimate;

% Get the number of outcomes (bins)
nOutcomes = qpParams.nOutcomes;

% Determine the number of bins to be reserved for upper and lower headroom
nLower = max([1 round(nOutcomes.*myQpfmriParams.headroom)]);
nUpper = max([1 round(nOutcomes.*myQpfmriParams.headroom)]);
nMid = nOutcomes - nLower - nUpper;

% Convert the adjusted BOLD amplitudes into scaled amplitudes (0-1)
scaledAmplitudes = adjustedAmplitudes./maxBOLDLatestGuess;

% Map the responses to binned outcomes
outcomes = 1+round(scaledAmplitudes.*nMid)+nLower;
outcomes(outcomes > nOutcomes)=nOutcomes;
outcomes(outcomes < 1)=1;

end % main function



