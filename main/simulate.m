function [psiParamsFit]=simulate(model_type, model_params, control_params, qpPres)
% model_type     -          string
%
% model_params   -          1xn vector of model parameters:
%                             DoE    (n=5): Sr, k1, k2, beta, sigma
%                             Watson (n=5): tau, kappa, zeta, beta, sigma
%
% control_params -         1x2 vector of control parameters:
%                             TRmsecs         - The TR of the BOLD fMRI measurement in 
%                                                   milliseconds. Needed in simulation mode to know how
%                                                   much to downsample the simulated signal. 
%                             trialLengthSecs - The number of seconds per trial
%
% qpPres         -         boolean:
%                           true  - run simulation with Q+ stimulus presentation
%                           false - run simulation without Q+ stimulus presentation
%
% Example: 
% model_type='doe';model_params=[0.9998 0.0132 0.7755 1 0];control_params=[800 12];qpPres='false';
% [psiParamsFit]=doeSimulate(model_params, control_params, qpPres);
%


if strcmp(model_type, 'doe')
    [psiParamsFit]=doeSimulate(num2str(model_params(1)), num2str(model_params(2)),...
    num2str(model_params(3)), num2str(model_params(4)), num2str(model_params(5)),...
    num2str(control_params(1)), num2str(control_params(2)), qpPres, '0');
elseif strcmp(model_type, 'wt')
    [psiParamsFit]=watsonSimulate(model_params, control_params, qpPres);
else
    disp('Invalid model')
end

end

