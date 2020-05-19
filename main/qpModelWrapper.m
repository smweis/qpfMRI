function predictedProportions = qpModelWrapper(model, stimulusValues, params, nOutcomesIn, headroomIn)
% Express a psychometric model as amplitude proportions
%
% Syntax:
%  predictedProportions = qpModelWrapper(model, stimulusValues, params, nOutcomesIn, headroomIn)
%
% Description:
%	This function maps the amplitude of the psychometric function response
%	to a discrete response within one of nCategory bins. The psychometric
%	model has amplitude values in the domain of [0 beta], where beta is the
%	parameter that controls maximum model amplitude. The params variable
%	defines the n params of the psychometric model, plus a final,
%	additional parameter that is the degree of Gaussian smoothing to be
%	applied to the response categories across the y-axis bins. The units of
%	this sigma value are nCategory bins of the response. So, a sigma of
%	0.25 indicates a Gaussian kernel that has a SD equal to 1/4 of the
%	respond range.
%
% Inputs:
%   model                 - Functional handle. The psychometric function to
%                           be mapped.
%   stimulusValues        - nx1 column vector of stimulus values for
%                           which the psychometric function ("model") will
%                           be evaluated and predicted proportions returned
%   params                - 1xn vector. The first n-1 values are the
%       	                parameters to be passed to the model. The nth
%                           value is the sigma of the Gaussian smoothing to
%                           apply across the category boundaries.
%   nOutcomesIn           - Scalar. Optional. Sets the number of bins into
%                           which the y-axis will be divided. Defaults to
%                           21 if not provided.
%   headroomIn            - Scalar. Optional. Determines the proportion
%                           of the nOutcomes to reserve above and below
%                           the minimum and maximum output of the Watson
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
    % Parameters of the Watson model

    % Define the psychometric function
    model = @watsonTemporalModel;

    % Define the model parameters
    tau = 1;	% time constant of the center filter (in msecs)
    kappa = 1.5;	% multiplier of the time-constant for the surround
	zeta = 1;	% multiplier that scales the amplitude of the surround
    beta = 1;   % maximum model output across frequencies

    % Plot the Watson TTF
    freqDomain = logspace(0,log10(100),100);
    figure
    subplot(2,1,1);
    semilogx(freqDomain,model(freqDomain,[tau kappa zeta beta]));
    xlabel('log Freq [Hz]');
    ylabel('Amplitude TTF [0-1]');
    title('Watson TTF');
    hold on

    % The frequency at which to obtain the predicted proportions
    freqHz = 40;

    % Indicate on the Watson TTF plot the frequency to model
    semilogx(freqHz,model(freqHz,[tau kappa zeta beta]),'*r');

    % Gaussian noise to be applied across the y-axis response categories
    sigma = 0.5;

    % Assemble the params with the noise modeled
    params = [tau kappa zeta beta sigma];

    % The number of bins into which to divide the response amplitude
    nOutcomes = 21;

    % Perform the calculation
    predictedProportions = qpModelWrapper(model,freqHz, params, nOutcomes);

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
if nargin >= 4
    nOutcomes = nOutcomesIn;
else
    nOutcomes = 21;
end

% Ensure that nOutcomes is odd
assert(mod(nOutcomes,2)==1);

% Set the headroom if undefined
if nargin >= 5
    headroom = headroomIn;
else
    headroom = 0.1;
end

% Determine the number of bins to be reserved for upper and lower headroom
nLower = round(nOutcomes.*headroom);
nUpper = round(nOutcomes.*headroom);
nMid = nOutcomes - nLower - nUpper;

% Obtain the model response values for the stimulus values
yVals = model(stimulusValues, params(1:end-1));

% Map the responses to categories
binAssignment = 1+round(yVals.*nMid)+nLower;
binAssignment(binAssignment > nOutcomes)=nOutcomes;

% Create a Gaussian kernel to reflect noise in the y-axis measurement
if params(end)==0
    gaussKernel = zeros(1,nOutcomes);
    gaussKernel(floor(nOutcomes/2)+1)=1;
else
    gaussKernel = gausswin(nOutcomes,nOutcomes/(10*params(end)))';
end

% Initialize a variable to hold the result
predictedProportions = zeros(length(stimulusValues),nOutcomes);

% Loop over the array of stimulusValues
for ii = 1:length(stimulusValues)
    
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


end % main function

