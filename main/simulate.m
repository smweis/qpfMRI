function []=simulate(model_type, model_params)
% model_type   -          string
%
% model_params -          1x5 vector of model parameters:
%                           Sr - ratio of the scale factors of the 
%                                exponentials
%                           k1, k2 - time constants of the exponentials
%                           beta - maximum amplitude of the response
%                           sigma - Standard deviation of the scaled (0-1) noise


    if strcmp(model_type, 'doe')
            validate_qpDoETFE_simulate(model_params);
    elseif strcmp(model_type, 'wt')
        validate_qpWatsonTFE_simulate(model_params);
    else
        disp('Invalid model')
    end

end