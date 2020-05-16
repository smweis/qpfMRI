function [paramNamesInOrder, qpPF, psiParamsDomainList] = checkModel(model, paramsDomain, varargin)
%checkModel will validate that the model and parameters specified match and
%are currently supported. 
% Inputs:
%   model                  - A function handle. This should be the
%                            'continuous function.' The Quest+ specific
%                            function will be defined in the model-specific
%                            code block. Currently supported:
%                              @doeTemporalModel
%                              @watsonTemporalModel
%   paramsDomain           - Struct consisting of upper bounds, lower
%                            bounds, and intervals for all necessary
%                            parameters. All models should have beta and
%                            sigma as parameters.
%                              DoE    (n=5): Sr, k1, k2, beta, sigma
%                              Watson (n=5): tau, kappa, zeta, beta, sigma
% Optional key/value pairs:
%   'headroom'              - Scalar: (Default = 0.1)
%                             The proportion of the nOutcomes from qpParams 
%                             that will be used as extra on top and bottom.
%	'nOutcomes'             - Integer (Default = 51)
%                             The number of outcome bins Q+ can assign a
%                             response to. The larger this is the slower
%                             Q+ will be. 
%Example
%{
model = @doeTemporalModel;

paramsDomain = struct;
paramsDomain.Sr = 0.899:0.025:1.099;
paramsDomain.k1 = 0.01:0.04:0.4;
paramsDomain.k2 = 0.01:0.04:0.4;
paramsDomain.beta = 0.8:0.1:1.4; % Amplitude of the scaled response; should converge to unity
paramsDomain.sigma = 0.3:0.2:1;	% Standard deviation of the scaled (0-1) noise

[paramNamesInOrder, qpPF, psiParamsDomainList] = checkModel(model,paramsDomain);

%}
%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));
p.addRequired('paramsDomain',@isstruct);

% Optional params
p.addParameter('headroom', 0.1, @isnumeric);
p.addParameter('nOutcomes',51,@isnumeric);

% Parse
p.parse( model, paramsDomain, varargin{:});

% P is used later
headroom = p.Results.headroom;
nOutcomes = p.Results.nOutcomes;

%% Model definitions
% To add a new model, copy the elseif block and add the new model name to 
% the models created variable, just below.
modelsCreated = {'doeTemporalModel','watsonTemporalModel'};

% A variable number of parameters are supported, but all models must
% contain a beta and a sigma parameter. 
s = functions(model);
assert(any(strcmp(modelsCreated,s.function)==1),'Model not defined');

% Model specific functions.
if contains(s.function,'doe')
    % This is important if the model expects input in a certain order, but
    % is required for current functionality.
    paramNamesInOrder = {'Sr', 'k1', 'k2', 'beta', 'sigma'};
    psiParamsDomainList = cell(1,length(paramNamesInOrder));
    for i = 1:length(paramNamesInOrder)
        % Check that all parameter names expected are present in the domain
        % and add the domain to myQpParams for Q+.
        assert(isfield(paramsDomain,paramNamesInOrder{i}),'Parameter missing or misnamed for DoE model.');
        psiParamsDomainList{i} = paramsDomain.(paramNamesInOrder{i});
    end
    % The QP form of the model should be specified below. 
    qpPF = @(f,p) qpDoETemporalModel(f,p,nOutcomes,headroom);
    
elseif contains(s.function,'watson')
    paramNamesInOrder = {'tau', 'kappa', 'zeta', 'beta', 'sigma'};
    for i = 1:length(paramNamesInOrder)
        assert(isfield(paramsDomain,paramNamesInOrder{i}),'Parameter missing or misnamed for Watson model.');
    end
    % Create an anonymous function from qpDoETemporalModel in which we
    % specify the number of outcomes for the y-axis response
    qpPF = @(f,p) qpWatsonTemporalModel(f,p,nOutcomes,headroom);
    
end

end
