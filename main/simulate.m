function [modelResponseStruct,thePacketOut,questData]=simulate(model_type, model_params, sim_type)
% model_type   -          string
%
% model_params -          1x5 vector of model parameters:
%                           Sr - ratio of the scale factors of the 
%                                exponentials
%                           k1, k2 - time constants of the exponentials
%                           beta - maximum amplitude of the response
%                           sigma - Standard deviation of the scaled (0-1) noise
%
% sim_type     -         boolean:
%                           true - run simulation without Q+ stimulus presentation
%                           false - run simulation with Q+ stimulus presentation


    if strcmp(model_type, 'doe')
        [modelResponseStruct,thePacketOut,questData]=validate_qpDoETFE_simulate(model_params, sim_type);
    elseif strcmp(model_type, 'wt')
        [modelResponseStruct,thePacketOut,questData]=validate_qpWatsonTFE_simulate(model_params);
    else
        disp('Invalid model')
    end

end


