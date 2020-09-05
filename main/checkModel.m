function [paramNamesInOrder, varargout] = checkModel(myQpfmriParams)
%checkModel will validate that the model and parameters specified match and
%are currently supported. 
% Inputs:
% Required Inputs
%   myQpfmriParams          - Struct. Set of parameters used for qpfmri.
%                             See qpfmriParams function for more details
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

myQpfmriParams = qpfmriParams(model,paramsDomain);

[paramNamesInOrder, qpPF, psiParamsDomainList] = checkModel(myQpfmriParams);

%}



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
s = functions(myQpfmriParams.model);
assert(any(strcmp(modelNames,s.function)==1),'Model not defined');

% Grab the right parameter names and qpPF function
paramNamesInOrder = modelParamNames.(s.function);

% Return the qpPF handle
varargout{1} = @(f,p) qpModelWrapper(myQpfmriParams.model,f,p,...,
    myQpfmriParams.nOutcomes,myQpfmriParams.headroom);

% Initialize the psiParamsDomain
varargout{2} = cell(1,length(paramNamesInOrder));

% Step through each parameter, making sure all fields are accounted
% for, then assign them to the cell.
for i = 1:length(paramNamesInOrder)
    assert(isfield(myQpfmriParams.paramsDomain,paramNamesInOrder{i}),'Parameter missing or misnamed for model.');
    varargout{2}{i} = myQpfmriParams.paramsDomain.(paramNamesInOrder{i});
end

end


