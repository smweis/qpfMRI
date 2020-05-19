function [paramNamesInOrder, varargout] = checkModel(model, varargin)
%checkModel will validate that the model and parameters specified match and
%are currently supported. 
% Inputs:
%   model                  - A function handle. This should be the
%                            'continuous function.' The Quest+ specific
%                            function will be defined in the model-specific
%                            code block. Currently supported:
%                              @doeTemporalModel
%                              @watsonTemporalModel
% Optional positional input:
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
% Outputs: 
%   paramNamesInOrder       - Cell. The list of parameter names in the order the
%                             model expects them to be. 
%
% Optional positional outputs:  
% 
%   qpPF                    - Cell. A cell array of the parameter domains
%   
%   psiParamsDomain         - function handle. The model for use with Q+.
%                             This can only be returned if paramsDomain
%                             is passed as an optional positional argument.
                           
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

% Just to get the paramNamesInOrder, can be called without paramsDomain
[paramNamesInOrder] = checkModel(model);

%}
%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));

% Optional positional argument, paramsDomain
p.addOptional('paramsDomain',{},@isstruct);

% Optional params
p.addParameter('headroom', 0.1, @isnumeric);
p.addParameter('nOutcomes',51,@isnumeric);

% Parse
p.parse( model, varargin{:});

% P is used later
paramsDomain = p.Results.paramsDomain;
headroom = p.Results.headroom;
nOutcomes = p.Results.nOutcomes;



%% ALL MODEL SPECIFIC STUFF HAS BEEN MOVED HERE! 
% To add a new model, copy the elseif block and add the new model name to 
% the models created variable, just below.

% Add all continuous model names to the list below. 
modelNames = {'doeTemporalModel','watsonTemporalModel','logistic'};

% Add all model parameters as a field in the struct below.
% This is important if the model expects input in a certain order, but
% is required for current functionality.
modelParamNames = struct;
modelParamNames.doeTemporalModel = {'Sr', 'k1', 'k2', 'beta', 'sigma'};
modelParamNames.watsonTemporalModel = {'tau', 'kappa', 'zeta', 'beta', 'sigma'};
modelParamNames.logistic = {'slope','semiSat','beta','sigma'};



%% Do the checks.

% A variable number of parameters are supported, but all models must
% contain a beta and a sigma parameter. 
s = functions(model);
assert(any(strcmp(modelNames,s.function)==1),'Model not defined');



% Grab the right parameter names and qpPF function
paramNamesInOrder = modelParamNames.(s.function);


% Return the qpPF handle
varargout{1} = @(f,p) qpModelWrapper(model,f,p,nOutcomes,headroom);



% If there's only one input argument but more than 2 output arguments, give
% a warning.
if nargin == 1
    createPsiParams = false;
    if nargout > 2
        warning('No psiParamsDomainList can be returned without psiParamDomains');
        varargout{2} = {};
    end
else
    createPsiParams = true;
end


if createPsiParams
    % Initialize the psiParamsDomain
    varargout{2} = cell(1,length(paramNamesInOrder));

    % Step through each parameter, making sure all fields are accounted
    % for, then assign them to the cell.
    for i = 1:length(paramNamesInOrder)
        assert(isfield(paramsDomain,paramNamesInOrder{i}),'Parameter missing or misnamed for Watson model.');
        varargout{2}{i} = paramsDomain.(paramNamesInOrder{i});
    end
    
end

end

