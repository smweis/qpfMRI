function y = nakaRushton(stimulusValue, params)
% Three-parameter Naka-Rushton function over the stimulus domain of 0-1
%
% Syntax:
%  y = nakaRushton(stimulusValue, params)
%
% Description:
%	It is valid over the domain of [0, 1]. The function is controlled by
%	three parameters, defined below.
%
%   The output domain of the model is [0, beta].
%
% Inputs:
%   stimulusValue         - 1xn vector that provides the stimulus values
%                           for which the model will be evaluated
%   params                - 1x2 vector of model parameters:
%                           	exponent
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
    y = nakaRushton(stimulusValues,params);
    plot(stimulusValues,y,'-k');
%}

% The valid stimulus domain for the model
stimulusDomain = [0 1];

% Sanity check the stimulus value input
if max(stimulusValue)>max(stimulusDomain) || min(stimulusValue)<min(stimulusDomain)
    error('The passed stimulus value is out of range for the model');
end

% Un-pack the passed parameters
exponent = params(1);
semiSat = params(2);
beta = params(3);

% Calculate the amplitudes
y = (stimulusValue.^exponent ./ (stimulusValue.^exponent + semiSat.^exponent));

% Scale to beta
y = beta .* y./max(y);

end % main function

