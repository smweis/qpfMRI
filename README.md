# qpFMRI
A MATLAB toolbox for adaptive stimulus presentation using [QUEST+](https://jov.arvojournals.org/article.aspx?articleid=2611972) in real-time fMRI experiments.  

## Usage
In MATLAB, run `tbUseProject('neurofeedback');` in the console and set your working directory to `toolboxes/qpfmri`. 
To begin, run `realtimequest(...)`. You will be prompted to select a global params file (`global.json`) and either generate or select a local params file (`qpParams.json`). 

This module reads in `fMRI_timeseries.txt` containing BOLD signal data and `actualStimuli.txt` containing presented stimuli. It then uses QUEST+ to dynamically suggest the next stimulus to present, and writes this to `suggestions.txt`. 

`realtimequest()` can be run in two modes: qp (1) or random (0); this is determined via the `qpPres` parameter. 

### Parameters
There are two ways to control the values of parameters in `qpfmri`: by file (default) or by command-line arguments. Changing parameter values by file is as simple as finding the desired parameter in `qpgetparams()` and editing its value. If you find yourself testing different parameter values often, setting the debug flag in `realtimequest.m` to `1` allows you to enter parameters and values as command line arguments (see example in `realtimequest.m`). 

### Important Locations
* `global.json`: `subjectProcessedPath/`
* `qpParams.json`: user-determined (`subjectProcessedPath/`)
* `fMRI_timeseries.txt`: `subjectProcessedPath/processed/run#/`
* `actualStimuli.txt`: `subjectProcessedPath/stims/run#/`
* `suggestions.txt`: `subjectProcessedPath/stims/run#/`
