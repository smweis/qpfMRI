function [model_params] = compiledSimulate(model_type, varargin)
%  formats CLI input for simulation
%
% model_type - string

p = inputParser;

if strcmp(model_type, 'doe')
    % Required input
    p.addRequired('model',@ischar);
    p.addRequired('Sr',@ischar);
    p.addRequired('k1',@ischar);
    p.addRequired('k2',@ischar);
    p.addRequired('beta',@ischar);
    p.addRequired('sigma',@ischar);
    p.addRequired('TR',@ischar);
    p.addRequired('trialLength',@ischar);
    p.addRequired('qpPres',@ischar);
    p.addRequired('outNum',@ischar);
    p.addRequired('seed',@ischar);

    % Replace any defaults then parse the input
    % Optional positional params
    p.addOptional('nTrials','30',@ischar);
    p.addOptional('stimulusStructDeltaT','100',@ischar);
    p.addOptional('maxBOLDSimulated','1.6',@ischar);
    p.addOptional('maxBOLD','1.0',@ischar);
    p.addOptional('baselineStimulus','0',@ischar);
    p.addOptional('maxBOLDStimulus','30',@ischar);
    p.addOptional('nOutcomes','51',@ischar);
    p.addOptional('headroom','.1',@ischar);
    
    p.parse(model_type, varargin{:});
    
    % build parameter structs
    model_params.Sr = p.Results.Sr;
    model_params.k1 = p.Results.k1;
    model_params.k2 = p.Results.k2;
    model_params.beta = p.Results.beta;
    model_params.sigma = p.Results.sigma;
    
    control_params.TR = p.Results.TR;
    control_params.trialLength = p.Results.trialLength;
    
end
simulate(p.Results.model, model_params, control_params, p.Results.qpPres);
end

