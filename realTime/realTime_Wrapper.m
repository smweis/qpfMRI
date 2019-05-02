

%% Explanation
% This script is a heavily modified version of runNeurofeedbackTemplate which will 
% grab and process real-time fMRI data at SC3T for use with Q+ and the TFE.


% Directory structure. Brackets = You create them. Parentheses = Script creates them.  
% [subject]
%   - [raw]             : T1w, scout EPI images, unprocessed
%   - [processed]       : T1w, scout EPI images after BET, as well as
%                         transformation matrices for EPI -> T1w -> MNI
%   
%    - (runX)           : Run-specific DICOMs. 


%% Either edit these variables here, or input them each time.

% subject - a string corresponding to the local directory on which your
%           subject-specific pre-processed data exists, and which the new data will
%           be stored. 
% run     - a string corresponding to the name of the run (can be an int)


subject = input('Subject number?','s');  
run = input('Which run?','s');
%sbrefQuestion = input('Is there an sbref (y or n)?','s');
%realOrTest = input('At SC3T (y or n)?','s');

realOrTest = 'n';


sbrefQuestion = 'y';

if mod(run,2) == 1
    ap_or_pa = 'PA';
else
    ap_or_pa = 'AP';
end



% set flags
if strcmp(sbrefQuestion,'y')
    registerToFirst = false;
    registerToSbref = true;
else
    registerToFirst = true;
    registerToSbref = false;
end



if strcmp(realOrTest,'y')
    atScanner = true;
else
    atScanner = false;
end

showFig = false;
checkForTrigger = true;



%% Get Relevant Paths

[subjectPath, scannerPathStem, codePath, scratchPath] = getPaths(subject);


% Set scanner path manually for this script:
scannerPath = strcat('/Users/nfuser/Documents/rtQuest/TOME_3021_rtMockScanner/simulatedScannerDirectory/run',run);


%{
% If we're at the scanner, get the most recently created folder on the scanner path.
if atScanner
    thisSessionPath = dir(scannerPathStem);
    thisSessionPathSorted = sortrows(struct2table(thisSessionPath),{'isdir','datenum'});
    scannerPath = strcat(table2cell(thisSessionPathSorted(end,'folder')), filesep, table2cell(thisSessionPathSorted(end,'name')));
else
    scannerPath = scannerPathStem;
end
%}


%% Register to First DICOM or SBREF

if registerToSbref
    %sbrefInput = input('Input the filename of the sbref including file type\n for example: 001_000013_000001.dcm','s');
    sbrefInput = strcat(ap_or_pa,'_run',run,'_SBRef.nii');
    sbref = [subjectPath filesep 'raw' filesep sbrefInput];
    
    [ap_or_pa,initialDirSize] = registerToFirstDicom(subject,subjectPath,run,scannerPath,codePath,sbref);
end

if registerToFirst
    [ap_or_pa,~] = registerToFirstDicom(subject,subjectPath,run,scannerPath,codePath);
end




%% Check for trigger

if checkForTrigger
    first_trigger_time = waitForTrigger;
end


%% Load the ROI

roiName = ['ROI_to_new',ap_or_pa,'_bin.nii.gz'];
roiPath = fullfile(subjectPath,'processed',strcat('run',run),roiName);
roiIndex = loadRoi(roiPath);


%% Main Neurofeedback Loop

% Initialize the main data struct;
mainData = struct;
mainData.acqTime = {}; % time at which the DICOM hit the local computer
mainData.dataTimepoint = {}; % time at which the DICOM was processed
mainData.dicomName = {}; % name of the DICOM
mainData.roiSignal = {}; % whatever signal is the output (default is mean)


% This script will check for a new DICOM, then call scripts that will
% convert it to NIFTI, and do some processing on the NIFTI.
% (Extract ROI, compute mean signal of the ROI).




i = 0;
j = 1;


while i < 10000000000
    i = i + 1;

    [mainData(j).acqTime,mainData(j).dataTimepoint,mainData(j).roiSignal,...
     initialDirSize, mainData(j).dicomName] = ...
     checkForNewDicom(scannerPath,roiIndex,initialDirSize,scratchPath);
    
    
    % Save the main data
    save(fullfile(subjectPath,'processed',strcat('run',run),strcat('mainDatarun',run)),'mainData');
 
 
    j = j + 1;

    pause(0.01);
    
end




