function [outcomes, modelResponseStruct, thePacket, adjustedAmplitudes] = tfeUpdate(thePacket, qpParams, stimulusVec, baselineStimulus, varargin)
% Returns the QP outcomes given a packet and a stimulus vector.
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
%   stimulusVec           - 1xk vector. Provides the numeric value for each
%                           stimulus.
%   baselineStimulus      - Scalar. Provides the numeric value for the
%                           baseline stimulus that is used for reference
%                           coding.
%
% Optional key/value pairs (used in fitting):
%  'headroom'             - Scalar. The proportion of the nOutcomes from 
%                           qpParams that will be used as extra on top and
%                           bottom.
%  'maxBOLD'              - Scalar. The value (in % change units) of the
%                           maximum expected response to a stimulus w.r.t.
%                           the response to the baseline stimulus.
%
% Optional key/value pairs (used in simulation):
%  'maxBOLDSimulated'     - Scalar. The value (in % change units) of the
%                           maximum expected response to a stimulus w.r.t.
%                           the response to the baseline stimulus.
%  'rngSeed'              - Struct. By passing a seed to the random
%                           number generator, calling function can ensure
%                           that this routine returns the same output for a
%                           given input in simulation model.
%  'noiseSD'              - Scalar. The amplitude of the noise added to the
%                           simulated BOLD fMRI signal in units of standard
%                           deviations.
%  'pinkNoise'            - Logical. If the noise should be modeled as
%                           pink (auto-correlated).
%  'TRmsecs'              - Scalar. The TR of the BOLD fMRI measurement in 
%                           msecs. Needed in simulation mode to know how
%                           much to downsample the simulated signal.
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
%
% Examples:
%{
    % SIMULATION MODE
    nTrials = 35;
    thePacket = createPacket('nTrials',nTrials);

    % Generate a random stimulus vector
    stimulusVec = randsample([0, 0, 1.875,3.75,7.5,10,15,20,30],nTrials,true);

    % Initialize some parameters to pass to tfeUpdate from Quest
    myQpParams = qpParams;
    % The number of outcome categories.
    myQpParams.nOutcomes = 51;
    
    % The headroom is the proportion of outcomes that are reserved above and
    % below the min and max output of the Watson model to account for noise
    headroom = .1;

    % Define binned psychometric fuction
    myQpParams.qpPF = @(f,p) qpWatsonTemporalModel(f,p,myQpParams.nOutcomes,headroom);

    % Create some simulatedPsiParams
    tau = 0.5:0.5:10;	% time constant of the center filter (in msecs)
    kappa = 0.5:0.25:3;	% multiplier of the time-constant for the surround
    zeta = 0:0.25:2;	% multiplier of the amplitude of the surround
    beta = 0.8:0.1:1;   % multiplier that maps watson 0-1 to BOLD % bins
    sigma = 0:0.25:2;	% width of the BOLD fMRI noise against the 0-1 y vals
    myQpParams.psiParamsDomainList = {tau, kappa, zeta, beta, sigma};
    simulatedPsiParams = [randsample(tau,1) randsample(kappa,1) randsample(zeta,1) randsample(beta,1) 1];

    % This is the continuous psychometric fuction
    myQpParams.continuousPF = @(f) watsonTemporalModel(f,simulatedPsiParams);

    % This is the veridical psychometric fuction in binned outcomes
    myQpParams.qpOutcomeF = @(f) qpSimulatedObserver(f,myQpParams.qpPF,simulatedPsiParams);
    
    % Identify which stimulus is the "baseline" stimulus
    baselineStimulus = 0;

    % Perform the simulation
    [binAssignment, modelResponseStruct, thePacketOut, adjustedAmplitudes] = tfeUpdate(thePacket, myQpParams, stimulusVec, baselineStimulus, 'headroom', headroom);

    % Plot the results
    figure
    subplot(2,1,1)
    plot(thePacketOut.response.timebase,thePacketOut.response.values,'.k');
    hold on
    plot(modelResponseStruct.timebase,modelResponseStruct.values,'-r');
    subplot(2,1,2)
    stimulusVecPlot = stimulusVec;
    stimulusVecPlot(stimulusVecPlot==0)=1;
    semilogx(stimulusVecPlot,binAssignment,'xk')
    refline(0,myQpParams.nOutcomes.*headroom);
    refline(0,myQpParams.nOutcomes-myQpParams.nOutcomes.*headroom);
    ylim([1 myQpParams.nOutcomes]);

%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('thePacket',@isstruct);
p.addRequired('qpParams',@isstruct);
p.addRequired('stimulusVec',@isnumeric);
p.addRequired('baselineStimulus',@isscalar);

% Optional params used in fitting
p.addParameter('headroom', 0.1, @isnumeric);
p.addParameter('maxBOLD', 1.0, @isscalar);

% Optional params used in simulation
p.addParameter('maxBOLDSimulated', 1.0, @isscalar);
p.addParameter('rngSeed',rng(1,'twister'),@isstruct);
p.addParameter('noiseSD',0.25, @isscalar);
p.addParameter('pinkNoise',1, @isnumeric);
p.addParameter('TRmsecs',800, @isnumeric);

% Parse
p.parse( thePacket, qpParams, stimulusVec, baselineStimulus, varargin{:});

% We need to have at least one "baseline" stimulus in the vector of
% stimuli to support reference-coding of the BOLD fMRI responses
if isempty(find(stimulusVec==p.Results.baselineStimulus, 1))
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


%% If response is empty, simulate it
if isempty(thePacket.response)
    
    % Initialize params0, which will allow us to create the forward model.
    params0 = tfeObj.defaultParams('defaultParamsInfo', defaultParamsInfo);
    params0.noiseSd = p.Results.noiseSD;
    params0.noiseInverseFrequencyPower = p.Results.pinkNoise;
    modelAmplitude = zeros(length(stimulusVec),1);
        
    % Obtain the continuous amplitude response
    for ii = 1:length(stimulusVec)        
        modelAmplitude(ii) = qpParams.continuousPF(stimulusVec(ii)) .* p.Results.maxBOLDSimulated;    
    end

    % We enforce reference coding, such that amplitude of response to the
    % stimuli is expressed relative to the response to the baseline
    % stimulus.
    modelAmplitude = modelAmplitude - mean(modelAmplitude==p.Results.baselineStimulus);

    % Place the responses in the paramMainMatrix
    params0.paramMainMatrix = modelAmplitude;
        
    % Lock the MATLAB random number generator to give us the same BOLD
    % noise on every iteration.
    rng(p.Results.rngSeed);
    
    % Create a simulated BOLD fMRI time series
    thePacket.response = tfeObj.computeResponse(params0,thePacket.stimulus,thePacket.kernel,'AddNoise',true);
    
    % Resample the response to the BOLD fMRI TR
    BOLDtimebase = min(thePacket.response.timebase):p.Results.TRmsecs:max(thePacket.response.timebase);
    thePacket.response= tfeObj.resampleTimebase(thePacket.response,BOLDtimebase);
    
    % Mean center the response
    thePacket.response.values = thePacket.response.values - mean(thePacket.response.values);
    
end


%% Fit the response
[params,~,modelResponseStruct] = tfeObj.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

% We engage in reference coding, such that the amplitude of any stimulus is
% expressed w.r.t. the "baseline" stimulus
adjustedAmplitudes = params.paramMainMatrix - ...
    mean(params.paramMainMatrix(stimulusVec==p.Results.baselineStimulus));

% Get the number of outcomes (bins)
nOutcomes = qpParams.nOutcomes;

% Determine the number of bins to be reserved for upper and lower headroom
nLower = round(nOutcomes.*p.Results.headroom);
nUpper = round(nOutcomes.*p.Results.headroom);
nMid = nOutcomes - nLower - nUpper;

% Convert the adjusted BOLD amplitudes into scaled amplitudes (0-1)
scaledAmplitudes = adjustedAmplitudes./p.Results.maxBOLD;

% Map the responses to binned outcomes
outcomes = 1+round(scaledAmplitudes.*nMid)+nLower;
outcomes(outcomes > nOutcomes)=nOutcomes;
outcomes(outcomes < 1)=1;

end % main function



