function qpfmriData = qpfmriParams(model,paramsDomain,varargin)
%qpfmriParams  Set user defined parameters for a QUEST+ fMRI run.
%
% Usage:
%     qpfmriData = qpfmriParams(varargin)
%
% Description:
%     Set the user defined parameters needed for a run of QUEST+ fMRI. These are
%     specified by key/value pairs.  The defaults are for a logistic regression example.
%
%     This works by allowing the user to pass a set of key value pairs for
%     each possible user defined parameter.
%
%
% Required inputs:
%   model                 - A function handle. This should be the
%                           'continuous function.' The Quest+ specific
%                           function will be defined in the model-specific
%                           code block. Currently supported:
%                             @doeTemporalModel
%                             @watsonTemporalModel
%                             @logistic
%   paramsDomain          - Struct consisting of upper bounds, lower
%                           bounds, and intervals for all necessary
%                           parameters. All models should have beta and
%                           sigma as parameters.
%                             DoE      (n=5): Sr, k1, k2, beta, sigma
%                             Watson   (n=5): tau, kappa, zeta, beta, sigma
%                             Logistic (n=4): slope, semiSat, beta, sigma
%
% Optional key/value pairs:
%
%   'qpPres'                -  Logical: (Default = false)
%                              true  - run simulation with Q+ stimulus choice
%                              false - run simulation with random stimulus
%                                      choice.
%   'simulatedPsiParams'    - Struct: (Default = randomly selected values
%                             from paramsDomain). The veridical parameters
%                             that are used for the forward model. Beta must
%                             be 1. 
%	'headroom'              - Scalar: (Default = 0.1)
%                             The proportion of the nOutcomes from qpParams 
%                             that will be used as extra on top and bottom.
%   'maxBOLDInitialGuess'   - Scalar: (Default = 1.0)
%                             The initial guess for the maximum BOLD value
%                             with respect to the baseline stimulus. 
%   'maxBOLDSimulated'      - Scalar: (Default = 1.5)
%                             The value (in % change units) of the
%                             maximum expected response to a stimulus w.r.t.
%                             the response to the baseline stimulus.
%	'seed'                  - No check (Default = 'choose')
%                             The value to initialize the rng. If 'choose'
%                             it will randomly initialize the seed using
%                             'shuffle'. 
%	'TR'                    - Integer: (Default = 800)
%                             Length of the time to repetition (TR) in
%                             milliseconds. 
%	'trialLength'           - Integer: (Default = 12)
%                             Length of one trial in seconds.
%	'outNum'                - String: (Default = 'test')
%                             Name of the output file (e.g., 'test.csv')
%	'outFolder'             - String: (Default = 'Results')
%                             Name of the output file (e.g., './Results')
%	'nTrials'               - Integer (Default = 10) 
%                             Number of trials to simulate. 
%	'stimulusStructDeltaT'  - Integer (Default = 100)
%                             Resolution of the stimulus struct in milliseconds
%                             (e.g., a value will be created every 100 ms).
%	'baselineStimulus'      - No check (Default = 0)
%                             tfeUpdate requires a baseline stimulus for
%                             which every value will be referenced. 
%	'maxBOLDStimulus'       - No check (Default = 15)
%                             It may help maxBOLD to require an early trial
%                             be a value we expect would result in a large
%                             BOLD response. 
%	'nOutcomes'             - Integer (Default = 51)
%                             The number of outcome bins Q+ can assign a
%                             response to. The larger this is the slower
%                             Q+ will be. 
%   'noiseSD'               - 1xn vector or scalar (Default = .1)
%                             Must be relative to maxBOLDSimulated. 
%                             In the absence of a simulatedPsiParam.sigma,
%                             noiseSD will allow the specification of a
%                             specific noise value (or selection from a
%                             vector of values). 
%	'stimulusDomain'        - Cell array (Default = a range of stimulus values).
%                             All possible stimulus values that can be
%                             assigned.
%   'stimulusDomainSpacing  - Default = 'lin'. 
%                             Whether the stimulusDomain is spaced linear
%                             or log.
% Outputs:
%     qpfmriData          Structure with one field each corresponding to the
%                        keys below. Each field has the same name as the
%                        key.

% 08/30/2020 smw Started on this.

%% Parse inputs and set defaults
%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));
p.addRequired('paramsDomain',@isstruct);

% Optional params
p.addParameter('qpPres',false,@islogical);
p.addParameter('simulatedPsiParams',{});
p.addParameter('headroom', 0.1, @isnumeric);
p.addParameter('maxBOLDInitialGuess', 1.0, @isscalar);
p.addParameter('maxBOLDSimulated', 1.5, @isscalar);
p.addParameter('TR',800, @isnumeric);
p.addParameter('trialLength',12, @isnumeric);
p.addParameter('outNum','test',@ischar);
p.addParameter('outFolder','Results',@ischar);
p.addParameter('seed','choose');
p.addParameter('nTrials',10,@isnumeric);
p.addParameter('stimulusStructDeltaT',100,@isnumeric);
p.addParameter('baselineStimulus','');
p.addParameter('maxBOLDStimulus','');
p.addParameter('nOutcomes',5,@isnumeric);
p.addParameter('noiseSD',.1,@isvector);
p.addParameter('stimulusDomain',{},@iscell);
p.addParameter('stimulusDomainSpacing','lin',@ischar);

% Parse
p.parse( model, paramsDomain, varargin{:});


%% Return structure
qpfmriData = p.Results;