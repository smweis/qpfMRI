function [myQpParams, questDataAtStart] = realTime_Init(modelType,modelParameters,varargin)
% Returns the QP outcomes given a packet and a stimulus vector.
%
% Syntax:
%  [myQpParams, questData] = realTime_Init(modelType,modelParameters,varargin)

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
%
% Optional key/value pairs (Q+ parameters):
%  'nOutcomes'            - Integer. Number of possible bins Q+ could use.
%  'stimulusDomain'       - Cell. 1xk vector of integers where k is the number of
%                           possible stimulus values that could be presented. 
%  'headroom'             - Scalar. The proportion of the nOutcomes from 
%                           qpParams that will be used as extra on top and
%                           bottom.
%
% Outputs:
%   myQpParams            - Struct. Parameter settings for use to initalize
%                           Q+. We may need some of these later.
%   questDataAtStart      - Initialized struct for use with Q+
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

[myQpParams, questDataAtStart] = realTime_qpInit(modelType,modelParameters);


% Using the Watson model
modelType = 'watson';
modelParameters = struct;
modelParameters.tau = 0.5:0.5:8;	% time constant of the center filter (in msecs)
modelParameters.kappa = 0.5:0.25:2;	% multiplier of the time-constant for the surround
modelParameters.zeta = 0:0.25:2;	% multiplier of the amplitude of the surround
modelParameters.beta = 0.5:0.2:2; % multiplier that maps watson 0-1 to BOLD % bins
modelParameters.sigma = 0:0.5:2;	% width of the BOLD fMRI noise against the 0-1 y vals

[myQpParams, questDataAtStart] = realTime_qpInit(modelType,modelParameters);

%}


%% Parse input
p = inputParser;

p.addRequired('modelType',@isstr);
p.addRequired('modelParameters',@isstruct);

% Optional params used for Q+
p.addParameter('nOutcomes', 51, @isnumeric);
p.addParameter('stimulusDomain', [0,1.875,3.75,7.5,15,30], @isvector);
p.addParameter('headroom', 0.1, @isnumeric);


% Parse
p.parse(modelType,modelParameters,varargin{:});





%% Set up Q+

% Get the default Q+ params
myQpParams = qpParams;

% Add the stimulus domain.
myQpParams.stimParamsDomainList = {p.Results.stimulusDomain};

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

