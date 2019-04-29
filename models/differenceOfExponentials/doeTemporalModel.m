function y = doeTemporalModel(frequenciesToModel, params)
% Difference of exponentials model
%
% Syntax:
%  y = doeTemporalModel(frequenciesToModel, params)
%
% Description:
%
%   Derrington, A. M., and P. Lennie. "Spatial and temporal contrast
%   sensitivities of neurones in lateral geniculate nucleus of macaque."
%   The Journal of physiology 357.1 (1984): 219-240.
%
% Inputs:
%   frequenciesToModel    - 1xn vector that provides the stimulus
%                           frequencies in Hz for which the model will be
%                           evaluated
%   params                - 1x4 vector of model parameters:
%                              Sr - ratio of the scale factors of the 
%                                   exponentials
%                          k1, k2 - time constants of the exponentials
%                            beta - maximum amplitude of the response
%
% Outputs:
%   y                     - 1xn vector of modeled amplitude values.
%
% Examples:
%{
    % Replicate Figure 7A
    freqHz = logspace(log10(0.1),log10(100),100);
    params = [517/513, 0.128, 0.135, 10];
    y = doeTemporalModel(freqHz,params);
    loglog(freqHz,y,'-k');
    hold on
    ylim([1 100]);
%}


% Define a frequency domain in Hz over which the model is defined. The
% maximum and minimum value of the y response should be contained within
% this range for any plausible set of parameters. We hard code a log-spaced
% range of 0 - 200 Hz here.
freqDomain = [0 logspace(log10(0.1),log10(200),100)];

% Sanity check the frequency input
if max(frequenciesToModel)>max(freqDomain) || min(frequenciesToModel)<min(freqDomain)
    error('The passed frequency is out of range for the model');
end

% Un-pack the passed parameters
Sr = params(1);
k1 = params(2);
k2 = params(3);
beta = params(4);

% Sanity check the params. If k1=k2, and Sr=1, then no output is possible
if and(k1==k2,Sr==1)
    error('The model is undefined when k1=k2, and Sr=1');
end

% Find the specific S values that produce the desired beta amplitude
F = @(w,S1,S2) (S1.*exp(-k1.*w)) - (S2.*exp(-k2.*w));
myObj = @(S1) beta - max(F(freqDomain,S1,S1./Sr));
S1 = fzero(myObj,1);
S2 = S1/Sr;

% Calculate the y values
y = F(frequenciesToModel,S1,S2);

end % main function
