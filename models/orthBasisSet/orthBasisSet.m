function y = orthBasisSet(stimulusValue, params)
% A three-parameter orthogonal basis set over the stimulus domain of 0-1
%
% Syntax:
%  y = orthBasisSet(stimulusValue, params)
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
%   params                - 1x4 vector of model parameters:
%                           	firstOrder
%                               secondOrder
%                               thirdOrder
%                               beta - max amplitude output
%
% Outputs:
%   y                     - 1xn vector of modeled amplitude values.
%
% Examples:
%{
    % Demonstrate basic output of the model
    stimulusValues = 0:0.1:1;
    colors = {'k','b','r','g','y'};
    figure;
    hold on;
    num=0;
    for a = -1:.2:1
        num = num+1;
        for b = -1:.2:1
            for c = 1
                params = [a b c 1];
                y = orthBasisSet(stimulusValues,params);
                plot(stimulusValues,y,'-');
            end
        end
    end
%}

% The valid stimulus domain for the model
stimulusDomain = [0 1];

% Sanity check the stimulus value input
if max(stimulusValue)>max(stimulusDomain) || min(stimulusValue)<min(stimulusDomain)
    error('The passed stimulus value is out of range for the model');
end

% Un-pack the passed parameters
firstOrder = params(1);
secondOrder = params(2);
thirdOrder = params(3);
beta = params(4);


% Transpose the stimulus values if they're a column rather than row.
if size(stimulusValue,1)>1
    stimulusValue = stimulusValue';
end

% Map the stimulus domain to -1:1
x = (stimulusValue.*2)-1;

% Polynomial function
function [ L2 ] = LegendreN( n, x )
% This function returns Legendre polynomial of degree n.
%
% Input parameters:
% 	n - order of polynomial
%	x - vector of time axis
%
% Returned value:
%	L2 - Legendre polynomial of degree n
L1 = zeros( 1, length(x) ); 
L2 = ones( 1, length(x) );
for i = 1:n
    L0 = L1;
    L1 = L2;
    L2 = ( (2.*i-1) .* x .* L1 - (i-1) .* L0) ./ i;
end
end

% First through third order polynomials
p1 = LegendreN(1,x);
p2 = LegendreN(2,x);
p3 = LegendreN(3,x);


% Calculate the amplitudes
y = firstOrder.*p1 + ...
    secondOrder.*p2 + ...
    thirdOrder.*p3;



% Normalize where 0 is the min of each function
stimDomainFine = stimulusDomain(1):.01:stimulusDomain(2);
xFine = (stimDomainFine.*2)-1;

p1Min = LegendreN(1,xFine);
p2Min = LegendreN(2,xFine);
p3Min = LegendreN(3,xFine);
 
fineY = firstOrder.*p1Min + ...
    secondOrder.*p2Min + ...
    thirdOrder.*p3Min;

% Normalize to be between 0 and 1
% Make sure we aren't dividing by zero (the max and min are the same).
if max(fineY) - min(fineY) > 1.e-6
    y = (y - min(fineY))./(max(fineY)-min(fineY));
end

% Scale to beta
y = beta .* y;

end % main function

