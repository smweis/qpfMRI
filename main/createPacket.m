function thePacket = createPacket(myQpfmriParams,varargin)
% Returns thePacket for use with tfeUpdate
%
% Syntax:
%  thePacket = createPacket();
%
% Description:
%	Generates a stimulusStruct and kernelStruct based on optional inputs
%	and prepares tfeObj and thePacket for use with tfeUpdate. 
%
% Inputs:
%   myQpfmriParams          - Struct. Set of parameters used for qpfmri.
%                             See qpfmriParams function for more details
%
% Optional key/value pairs:
%   'verbose'               - How talkative. 
%                             Default - False
% Outputs:
%   thePacket            - Structure.
%
% Examples:
%{
   thePacket = createPacket(myQpfmriParams)
%}


%% Parse input
p = inputParser;

% Required input
p.addRequired('myQpfmriParams',@isstruct);
% Optional params
p.addParameter('verbose', false, @islogical);

% Parse and check the parameters
p.parse( myQpfmriParams, varargin{:});


%% Temporal domain of the stimulus
deltaT = myQpfmriParams.stimulusStructDeltaT; % in msecs
totalTime = myQpfmriParams.nTrials*myQpfmriParams.trialLength*1000; % in msecs.
eventDuration = myQpfmriParams.trialLength*1000; % block duration in msecs

% Define the timebase
stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
nTimeSamples = size(stimulusStruct.timebase,2);


% This loop will create a stimulus struct that has the property that each 
% non-baseline trial will create its own regressor (row). Each baseline
% trial will be added to the first regressor.
for ii=1:(myQpfmriParams.nTrials)
    stimulusStruct.values(ii,:)=zeros(1,nTimeSamples);
    stimulusStruct.values(ii,(ii-1)*eventDuration/deltaT+1:ii*eventDuration/deltaT)=1;
end

%% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets

kernelStruct.timebase=linspace(0,15999,16000);

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
    gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;

% Normalize the kernel to have unit amplitude
[ kernelStruct ] = normalizeKernelArea( kernelStruct );

%% Construct a packet and model params
thePacket.stimulus = stimulusStruct;
thePacket.response = [];
thePacket.kernel = kernelStruct;
thePacket.metaData = [];


end