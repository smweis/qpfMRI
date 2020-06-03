function [domain] = makeDomain(lower,upper,nDivisions,varargin)
%[domain] = makeDomain(lower,upper,nDivisions,varargin)
%
% Create the parameter domain
% 
%
% Optional key/value pairs
%   'spacing' - String. Default = 'lin'
%               What spacing to use for the variable. Can be: 
%               'lin'  - linear spacing
%               'log'  - log spacing
%               'zeno' - asymptotically approach a halfway between lower /
%                        upper.
% 

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; 
p.addRequired('lower',@isscalar);
p.addRequired('upper',@isscalar);
p.addRequired('nDivisions',@isnumeric);

p.addParameter('spacing','lin',@ischar);

p.parse(lower,upper,nDivisions,varargin{:});



%% Do checks

assert(lower<upper,'Improper lower and upper bounds specified.');

if contains(p.Results.spacing,'lin')
    domain = linspace(lower,upper,nDivisions);
elseif contains(p.Results.spacing,'log')
    domain = logspace(lower,upper,nDivisions);
elseif contains(p.Results.spacing,'zeno')
    domain = zeros(1,nDivisions);
    midpoint = ((upper - lower)./2) + lower;
    assert(mod(nDivisions,2)==1,'Spacing specified as zeno, but nDivisions must be odd.');
    
    % This will assign the zeno spaced parameters
    for i = 1:(nDivisions-1)/2
        if i == 1
            thisLowerValue = lower;
            thisUpperValue = upper;
        else
            thisLowerValue = thisLowerValue + ((midpoint - lower)/2^(i-1));
            thisUpperValue = thisUpperValue + ((midpoint - upper)/2^(i-1));
        end
        domain(i) = thisLowerValue;
        domain(end + 1 - i) = thisUpperValue;
    end
    domain((nDivisions+1)/2) = midpoint;
else
    error('Improper spacing key %s. Log or linear spacing supported',p.Results.spacing);
end

 
