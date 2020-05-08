function [plotParams] = plotParamsDomain(model, params)
% model_type     -          functionHandle.
%
% model_params   -          1xn vector of model parameters:
%                             DoE    (n=5): Sr, k1, k2, beta, sigma
%                             Watson (n=5): tau, kappa, zeta, beta, sigma
%
% Example: 
%{
% Params domain to plot
params = struct;
params.Sr = 0.899:0.025:1.099;
params.k1 = 0.01:0.04:0.4;
params.k2 = 0.01:0.04:0.4;
params.beta = 1;
params.sigma = 0;

model = @doeTemporalModel;

plotParamsDomain(model, params);
%
%}

% Model specific processing. This should be functionalized. 
modelsCreated = {'doeTemporalModel','watsonTemporalModel'};
s = functions(model);
assert(any(strcmp(modelsCreated,s.function)==1),'Model not defined');

if contains(s.function,'doe')
    assert(isfield(params,'Sr'),'Sr missing or misnamed for DoE model.');
    assert(isfield(params,'k1'),'k1 missing or misnamed for DoE model.');
    assert(isfield(params,'k2'),'k2 missing or misnamed for DoE model.');
elseif contains(s.function,'watson')
    assert(isfield(params,'tau'),'tau missing or misnamed for Watson model.');
    assert(isfield(params,'kappa'),'kappa missing or misnamed for Watson model.');
    assert(isfield(params,'zeta'),'zeta missing or misnamed for Watson model.');
end

stimulusFreqHzFine = logspace(log10(0.01),log10(100),100);

parameters = fieldnames(params);

paramSpaceSize = zeros(1,length(parameters));

for p = 1:length(parameters)
    paramSpaceSize(p) = length(params.(parameters{p}));
end

figure

plotParams = zeros(length(parameters));
for p = 1:length(parameters)
    for i = 1:length(params.(parameters{p}))
        plotParams(p) = params.(parameters{p})(i);
        semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,params),'-k');
    end
end

model(stimulusFreqHzFine,[params.Sr(1) params.k1(1) params.k2(1) params.beta(1)])
semilogx(stimulusFreqHzFine,doeTemporalModel(stimulusFreqHzFine,params),'-k');
hold on
