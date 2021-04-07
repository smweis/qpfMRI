function [myQpfmriParams, varargout] = checkModel(myQpfmriParams)
%% [myQpfmriParams, varargout] = checkModel(myQpfmriParams)
% Helper function to ensure that the model and parameters specified match and
% are currently supported. 
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
myQpfmriParams = qpfmriParams(model,paramsDomain);

[myQpfmriParams, myQpParams.qpPF, myQpParams.psiParamsDomainList] = checkModel(myQpfmriParams);

%}



%% ALL MODEL SPECIFIC STUFF HAS BEEN MOVED HERE! 
% To add a new model, copy the elseif block and add the new model name to 
% the models created variable, just below.

% Add all continuous model names to the list below. 
modelNames = {'doeTemporalModel','watsonTemporalModel','logistic','nakaRushton'};

% Add all model parameters as a field in the struct below.
% This is important if the model expects input in a certain order, but
% is required for current functionality.
modelParamNames = struct;
modelParamNames.doeTemporalModel = {'Sr', 'k1', 'k2', 'beta', 'sigma'};
modelParamNames.watsonTemporalModel = {'tau', 'kappa', 'zeta', 'beta', 'sigma'};
modelParamNames.logistic = {'slope','semiSat','beta','sigma'};
modelParamNames.nakaRushton = {'exponent','semiSat','beta','sigma'};


%% Do the checks.

% A variable number of parameters are supported, but all models must
% contain a beta and a sigma parameter. 
s = functions(myQpfmriParams.model);
assert(any(strcmp(modelNames,s.function)==1),'Model not defined');

% Grab the right parameter names and qpPF function
myQpfmriParams.paramNamesInOrder = modelParamNames.(s.function);

% Add betaIndex and sigmaIndex to myQpfmriParams
myQpfmriParams.betaIndex = find(strcmp(myQpfmriParams.paramNamesInOrder,'beta'));
myQpfmriParams.sigmaIndex = find(strcmp(myQpfmriParams.paramNamesInOrder,'sigma'));

% Return the qpPF handle
varargout{1} = @(f,p) qpModelWrapper(myQpfmriParams.model,f,p,...,
    myQpfmriParams.nOutcomes,myQpfmriParams.headroom,myQpfmriParams.betaIndex,myQpfmriParams.sigmaIndex);

% Initialize the psiParamsDomain
varargout{2} = cell(1,length(myQpfmriParams.paramNamesInOrder));

% Step through each parameter, making sure all fields are accounted
% for, then assign them to the cell.
for i = 1:length(myQpfmriParams.paramNamesInOrder)
    assert(isfield(myQpfmriParams.paramsDomain,myQpfmriParams.paramNamesInOrder{i}),'Parameter missing or misnamed for model.');
    varargout{2}{i} = myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{i});
end

end


