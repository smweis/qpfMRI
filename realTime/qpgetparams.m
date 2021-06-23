function [subject,myQpfmriParams,myQpParams] = qpgetparams()
% Sets up moduel parameters

% Load global params
disp('Select global params file');
[file,path] = uigetfile('.json');
fid = fopen(fullfile(path,file),'rt');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
globalParams = jsondecode(strrep(str,'\','\\')); % add escape chars

subject = globalParams.subject;
run = globalParams.run;

loadParams = input('Choose params file [1] or generate one [2]: ','s');
if strcmp(loadParams,'1')
    % Load  parameters
    [file,path] = uigetfile('.json');
    fid = fopen(fullfile(path,file),'rt');
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    p = jsondecode(str);
    disp(strcat('Loaded ',file));

    p.stimulusDomain = {p.stimulusDomain};
   
elseif strcmp(loadParams,'2')
    % Generate params
    
    % Provide a model handle
    p.model = @logistic;
    
    % Specify the parameter domain. Each value must correspond to a parameter
    % expected by the mymodel.
    
    p.paramsDomain = struct;
    p.paramsDomain.slope = makeDomain(-1.2,-.2,10,'spacing','log');
    p.paramsDomain.semiSat = makeDomain(.01,1,10);
    p.paramsDomain.beta = makeDomain(.75,1.25,11,'spacing','zeno');
    
    % Sigma in the parameter domain is searching for noiseSD
    p.paramsDomain.sigma = makeDomain(0.01,1,10);
    
    % Specify a stimulus domain and whether it spaced linear or log.
    p.stimulusDomain = {makeDomain(.5,1,12)};
    p.stimulusDomainSpacing = 'log';
    
    % Number of trials to run and TR.
    p.nTrials = 19;
    p.TRmsecs = 800;
    
    % Allow Q+ to control the stimuli or not (false).
    p.qpPres = true;
    
    % Set the number of outcome categories / bins.
    p.nOutcomes = 15;
    
    % Do you want to see plots?
    %showPlots = true;
    
    %How long the trials are (in seconds).
    p.trialLength = 12;
    
    p.outNum = 0;
    
    p.model = func2str(p.model);
    
    disp('Choose a location for the params file [subjectProcessedPath/]');
    selPath = uigetdir;
    filename = strcat('qpParams_',subject,run,'.json');
    fid = fopen(fullfile(selPath,filename'),'wt');
    fprintf(fid,jsonencode(p));
    fclose(fid);
else
    error('Invalid choice');
end
% formatting
p.model = str2func(p.model);
[myQpfmriParams,myQpParams] = qpfmriParams(p.model,p.paramsDomain,'qpPres',p.qpPres,...,
    'stimulusDomain',p.stimulusDomain,'stimulusDomainSpacing',p.stimulusDomainSpacing,...,
    'nTrials',p.nTrials,'trialLength',p.trialLength,'nOutcomes',p.nOutcomes,'TR',p.TRmsecs,'outNum',p.outNum);

myQpfmriParams.outNum = str2double(run);




end


