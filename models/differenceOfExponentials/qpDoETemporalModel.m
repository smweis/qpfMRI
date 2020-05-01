function predictedProportions = qpDoETemporalModel(frequenciesToModel, params, nOutcomesIn, headroomIn)
% Express the DoE model TTF as amplitude proportions
%
% Syntax:
%  predictedProportions = qpDoETemporalModel(frequenciesToModel, params, nOutcomesIn, headroomIn)
%
% Description:
%	This function maps the 0-1 amplitude of the DoE TTF response to a
%	discrete response within one of nCategory bins. The shape of the DoE
%	TTF is determined by the first three elements of the params variable.
%	The fourth element of params is the degree of Gaussian smoothing to be
%	applied to the response categories across the y-axis. The units of this
%	sigma value are nCategory bins of the response. So, a sigma of 0.25
%	indicates a Gaussian kernel that has a SD equal to 1/4 of the respond
%	range.
%
% Inputs:
%   frequenciesToModel    - nx1 column vector of frequencies (in Hz) for 
%                           which the DoE TTF will be evaluated and 
%                           predicted proportions returned
%   params                - 1x5 vector. The first three values correspond
%                           to the three parameters of doeTemporalModel.
%                           The second to last parameter is a scaling
%                           factor, mapping percent BOLD to a 0-1 range.
%                           The last parameter is the sigma of the Gaussian
%                           smoothing to apply across the category
%                           boundaries.
%   nOutcomesIn           - Scalar. Optional. Sets the number of bins into
%                           which the y-axis will be divided. Defaults to
%                           21 if not provided.
%   headroomIn            - Scalar. Optional. Determines the proportion
%                           of the nOutcomes to reserve above and below
%                           the minimum and maximum output of the DoE
%                           model. Defaults to 0.1, which means that 20% of
%                           the nOutcomes range will correspond to
%                           response values that are less than zero or
%                           greater than 1.
% Outputs:
%   predictedProportions  - An nFrequencies x nOutcomes matrix that
%                           provides for each of the frequencies to model
%                           the probability that a measured response will
%                           fall in a given amplitude bin.
%
% Examples:
%{
    % Demonstrate the conversion of a single modeled frequency into a
    % predicted proportion vector under the control of smoothing noise
    % Parameters of the DoE model
    Sr = 1;
    k1 = 0.128;	% multiplier of the time-constant for the surround
    k2 = 0.135;
    beta = 1;   % maximum model output across frequencies

    % Plot the doe TTF
    freqDomain = logspace(0,log10(100),100);
    figure
    subplot(2,1,1);
    semilogx(freqDomain,doeTemporalModel(freqDomain,[Sr k1 k2 beta]));
    xlabel('log Freq [Hz]');
    ylabel('Amplitude TTF [0-1]');
    title('DoE TTF');
    hold on

    % The frequency at which to obtain the predicted proportions
    freqHz = 20;

    % Indicate on the DoE TTF plot the frequency to model
    semilogx(freqHz,doeTemporalModel(freqHz,[Sr k1 k2 beta]),'*r');

    % Gaussian noise to be applied across the y-axis response categories
    sigma = 0.5;

    % Assemble the params with the noise modeled
    params = [Sr k1 k2 beta sigma];

    % The number of bins into which to divide the response amplitude
    nOutcomes = 21;

    % Perform the calculation
    predictedProportions = qpDoETemporalModel(freqHz, params, nOutcomes);

    % Plot the predicted proportions
    subplot(2,1,2);
    plot(predictedProportions,1:nOutcomes,'*k')
    xlim([0 1]);
    ylim([1 nOutcomes]);
    xlabel('predicted proportion [0-1]');
    ylabel('Amplitude response bin');
    title('Predicted proportions');
%}

% The number of bins into which we will divide the range of y-axis response
% values. Either passed or set as a default

if nargin >= 3
    nOutcomes = nOutcomesIn;
else
    nOutcomes = 21;
end

% Ensure that nOutcomes is odd
assert(mod(nOutcomes,2)==1);

% Set the headroom if undefined
if nargin >= 4
    headroom = headroomIn;
else
    headroom = 0.1;
end

%% Params to vars
Sr = params(1);
k1 = params(2);
k2 = params(3);
beta = params(4);
sigma = params(5);	% width of the BOLD fMRI noise against the 0-1 y vals

% Determine the number of bins to be reserved for upper and lower headroom
nLower = round(nOutcomes.*headroom);
nUpper = round(nOutcomes.*headroom);
nMid = nOutcomes - nLower - nUpper;

% Obtain the response values for the frequencies to be modeled
yVals = doeTemporalModel(frequenciesToModel, [Sr k1 k2 beta]);

% Map the responses to categories
binAssignment = 1+round(yVals.*nMid)+nLower;
binAssignment(binAssignment > nOutcomes)=nOutcomes;
binAssignment(binAssignment < 1)=1;

% Create a Gaussian kernel to reflect noise in the y-axis measurement
if params(end)==0
    gaussKernel = zeros(1,nOutcomes);
    gaussKernel(floor(nOutcomes/2)+1)=1;
else
    gaussKernel = gausswin(nOutcomes,nOutcomes/(10*sigma))';
end

% Initialize a variable to hold the result
predictedProportions = zeros(length(frequenciesToModel),nOutcomes);

% Loop over the array of frequenciesToModel
for ii = 1:length(frequenciesToModel)

    % Assign the bin
    predictedProportions(ii,binAssignment(ii))=1;

    % Apply smoothing across the proportions of response predicted for this
    % frequency. This models the effect of noise in the measurement of
    % response amplitude
    predictedProportions(ii,:) = conv(predictedProportions(ii,:),gaussKernel,'same');
    
    % Ensure that the sum of probabilities across categories for a given
    % frequency is unity
    predictedProportions(ii,:) = predictedProportions(ii,:)/sum(predictedProportions(ii,:));
        
end % loop over frequencies to model

% log likelihood calculations upon these proportions are unhappy
% with zeros. We set them here instead to realmin
predictedProportions(predictedProportions==0)=realmin;

end % main function

