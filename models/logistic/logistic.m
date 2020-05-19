function y = logistic(stimulusValue, params)
% Simple two-parameter logisitic function over the stimulus domain of 0-1
%
% Syntax:
%  y = logistic(stimulusValue, params)
%
% Description:
%	A two-parameter logisitic model. It is valid over the domain of [0, 1].
%	stimulus value. The slope parameter is defined over [-1, 1], and the
%	semi-saturation parameter between [0, 1].
%
%   The output domain of the model is [0, beta].
%
% Inputs:
%   stimulusValue         - 1xn vector that provides the stimulus values
%                           for which the model will be evaluated
%   params                - 1x2 vector of model parameters:
%                           	slope
%                               semiSat
%                               beta - max amplitude output
%
% Outputs:
%   y                     - 1xn vector of modeled amplitude values.
%
% Examples:
%{
    % Demonstrate basic output of the model
    stimulusValues = 0:0.1:1;
    params = [0.2 0.5 1];
    y = logistic(stimulusValues,params);
    plot(stimulusValues,y,'-k');
%}

% The valid stimulus domain for the model
stimulusDomain = [0 1];

% Sanity check the stimulus value input
if max(stimulusValue)>max(stimulusDomain) || min(stimulusValue)<min(stimulusDomain)
    error('The passed stimulus value is out of range for the model');
end

% Un-pack the passed parameters
slope = params(1);
semiSat = params(2);
beta = params(3);

% Calculate the amplitudes
y = beta.*(1-1./(1+exp(slope.*(100.*(stimulusValue-semiSat)))));

end % main function

