function [myQpfmriParams] = createRandomModelParams(myQpfmriParams)
%qpfmriParams  Simulate a set of random model parameters from the parameter domains.
%
% Usage:
%     [simulatedPsiParams] = createRandomModelParams(myQpfmriParams)
%
% Description:
%     Randomly generate a set of model parameters based on the min and max bounds of
%     the parameter domains.
%
%
% Required inputs:
%   myQpfmriParams        - Struct. See qpfmriParams
%
% Optional key/value pairs:
%
% Outputs:
%   myQpfmriParams      - Struct with simulatedPsiParams generated. See qpfmriParams
%
% 04/22/2021 smw Started on this.
%Examples: 
%{
% Provide a model handle
model = @nakaRushton;

paramsDomain = struct;
paramsDomain.exponent = makeDomain(-1.2,-.2,10,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.beta = makeDomain(.75,1.25,11,'spacing','zeno');

paramsDomain.sigma = makeDomain(.5,4,8);

% createRandomModelParams requires generating qpfmriParams
[myQpfmriParams,myQpParams] = qpfmriParams(model,paramsDomain);

% Print the simulatedPsiParams
myQpfmriParams.simulatedPsiParams

%}

%% Parse inputs and set defaults
%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('myQpfmriParams',@isstruct);

% Parse
p.parse( myQpfmriParams );

myQpfmriParams.simulatedPsiParams = zeros(1,length(myQpfmriParams.paramNamesInOrder));
stillSearching = true;
while stillSearching
    for i = 1:length(myQpfmriParams.paramNamesInOrder)
        % Select a random value for between the minimum parameter value
        % and the maximum parameter value. 
        minBound = min(myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{i}));
        maxBound = max(myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{i}));
        myQpfmriParams.simulatedPsiParams(i) = minBound + (maxBound-minBound).*rand(1,1);
    end

    % Beta is always one for simulations
    myQpfmriParams.simulatedPsiParams(myQpfmriParams.betaIndex) = 1;

    % Simulated noise is selected from a random sample of noiseSD
    myQpfmriParams.simulatedPsiParams(myQpfmriParams.sigmaIndex) = randsample(myQpfmriParams.noiseSD,1);

    if strcmp(func2str(myQpfmriParams.model),'logistic')
        if abs(myQpfmriParams.model(myQpfmriParams.baselineStimulus,myQpfmriParams.simulatedPsiParams)) < myQpfmriParams.simulatedPsiParams(myQpfmriParams.betaIndex)/10000 && ...
                abs(myQpfmriParams.model(myQpfmriParams.maxBOLDStimulus,myQpfmriParams.simulatedPsiParams)) < 1 && ...
                abs(myQpfmriParams.model(myQpfmriParams.maxBOLDStimulus,myQpfmriParams.simulatedPsiParams)) > .99
            stillSearching = false;
        end
    else
        stillSearching = false;
    end
end

% Create baseline stimulus and maxBOLD stimulus
[~, baselineStimulusIndex] = min(myQpfmriParams.model(myQpfmriParams.stimulusDomain{:},myQpfmriParams.simulatedPsiParams));
myQpfmriParams.baselineStimulus = myQpfmriParams.stimulusDomain{:}(baselineStimulusIndex);
[~, maxBOLDStimulusIndex] = max(myQpfmriParams.model(myQpfmriParams.stimulusDomain{:},myQpfmriParams.simulatedPsiParams));
myQpfmriParams.maxBOLDStimulus = myQpfmriParams.stimulusDomain{:}(maxBOLDStimulusIndex);

end

