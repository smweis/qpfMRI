function summarizeAndPlotSimulations(model,sameParams,stimDomain,baseline,reloadData,dirName,factorName,paramNamesToPlot,paramDomains,nParamsForHist)
%% summarizeAndPlotSimulations(model,sameParams,stimDomain,baseline,reloadData,dirName,factorName,paramNamesToPlot,paramDomains,nParamsForHist)

% Plot and print some summary statistics for simulations results


% Syntax:
%  [psiParamsFit]=simulate(model, paramsDomain, varargin)
%
% Description:
%	Takes in a model and a possible set of parameters and whether or not Q+ 
%   is in control of things flag. 
%
% Inputs:
%   model                 - A function handle. This should be the
%                           'continuous function.' The Quest+ specific
%                           function will be defined in the model-specific
%                           code block. Currently supported:
%                             @doeTemporalModel
%                             @watsonTemporalModel
%   paramsDomain          - Struct consisting of upper bounds, lower
%                           bounds, and intervals for all necessary
%                           parameters. All models should have beta and
%                           sigma as parameters.
%                             DoE    (n=5): Sr, k1, k2, beta, sigma
%                             Watson (n=5): tau, kappa, zeta, beta, sigma
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
%   'maxBOLD'               - Scalar: (Default = 1.0)
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
%	'questDataCopy'         - Struct (Default is empty)
%                             If it is not necessary to reinitialize Q+
%                             (questData untrained is in memory), you can
%                             pass the initialized questDataCopy to
%                             simulate to save the time initializing. Note,
%                             any change to paramsDomain requires
%                             re-initialization. Be cautious changes other
%                             options without re-initializing. 
% Optional key/value pairs (used in plotting):
%   'showPlots'             - Logical: (Default = false)
%                             Whether to show plots.
%   'figWidth'              - Integer (Default = 900)
%                             Width of figure window size.
%   'figHeight'             - Integer (Default = 900)
%                             Height of figure window size.

% Outputs:
%   psiParamsFit          - 1xn vector returning the BADS best fit for the
%                           parameters
%   maxBOLD               - Scalar. Best estimate at the maximum BOLD
%                           value.
%   questDataCopy         - Struct. Copy of initialized questData.

%Example: 
%{
---------------------------------------------------------------------------
% Example 1:

model = @logistic;
sameParams = true; % Are the veridical model params the same for all parameters?
stimDomain = linspace(.01,1,25); % What's the stimDomain?
baseline = .01;
reloadData = true; % Do we need to reload the data into memory?
dirStem = pwd;
dirName = fullfile(dirStem,'logisticResultsParamSet1');

% What factor do you want to look at? 
factorName = 'sigmaSim';
%factorName = 'nOutcomes';

% How many parameters to look at in the main figure?
paramNamesToPlot = {'slope','semiSat','maxBOLD','sigma'};
paramDomains = {linspace(0,1,100),linspace(0,1,100),linspace(0,2.5,100),linspace(.1,9,100)};
nParamsForHist = length(paramNamesToPlot);

% Example 2: Make figures for all parameter sets: 


model = @logistic;
sameParams = true; % Are the veridical model params the same for all parameters?
stimDomain = linspace(.01,1,25); % What's the stimDomain?
baseline = .01;
reloadData = true; % Do we need to reload the data into memory?

% What factor do you want to look at? 
factorName = 'sigmaSim';
%factorName = 'nOutcomes';

% How many parameters to look at in the main figure?
paramNamesToPlot = {'slope','semiSat','maxBOLD','sigma'};
paramDomains = {linspace(0,1,100),linspace(0,1,100),linspace(0,2.5,100),linspace(.1,9,100)};
nParamsForHist = length(paramNamesToPlot);


a = dir('logistic*');

for name = 1:length(a) 
    summarizeAndPlotSimulations(model,sameParams,stimDomain,baseline,reloadData,...,
    a(name).name,factorName,paramNamesToPlot,paramDomains,nParamsForHist)
end
%}

%% Check model, grab parameter names, load and sanity check data
paramNamesInOrder = checkModel(model);

% If 'data' not in memory, run the function to extract the data from the directory
if reloadData
    fullDirName = fullfile(dirName,'Results');
    [data] = loadSimulatedData(fullDirName,model,[1,2]);
end

%% Print averages by Q+ selection and a factor of your choice

% If params are all the same, we can use this: 
if sameParams
    sampleSimulatedParams = struct;
    sampleSimulatedParamsVec = zeros(1,length(paramNamesInOrder));
    fprintf('One veridical parameter set is used:\n');
    for j = 1:nParamsForHist
        sampleSimulatedParams.(paramNamesToPlot{j}) = data.(strcat(paramNamesToPlot{j},'Sim'))(1);
        fprintf('Simulated parameter "%s" = %.03f\n',paramNamesToPlot{j},data.(strcat(paramNamesToPlot{j},'Sim'))(1));
    end
    for j = 1:length(paramNamesInOrder)
        sampleSimulatedParamsVec(j) = data.(strcat(paramNamesInOrder{j},'Sim'))(1);
        sampleMaxBOLD = data.maxBOLDSim(1);
    end
    predictedRelativeResponse = makePredicted(model, stimDomain, sampleSimulatedParamsVec, baseline, sampleMaxBOLD);
else
    error('multiple params not handled at this point');
end

% Print the average result for Q+/random and the given factor. 
avgResults = varfun(@mean,data,'InputVariables',[paramNamesInOrder 'maxBOLD'],...
       'GroupingVariables',{factorName,'qpPres'})



%% Plotting colors
% Enter Colors for plot
lightRandColor = '#3D92C9';
darkRandColor = '#165172';
lightQPColor = '#F99C16';
darkQPColor = '#FA4515';
veridicalColor = '#000000';



% Detect number of levels for the given factor name:
[~, ind] = unique(data.(factorName), 'rows');
qpAverageParams = zeros(length(ind),length(paramNamesInOrder));
randomAverageParams = zeros(length(ind),length(paramNamesInOrder));

% Finding the relation between sigma, maxBOLD, and nOutcomes
%data.sbnConstant = (data.maxBOLD ./ data.sigma);


qPlusPanels = zeros(nParamsForHist,1);
randPanels = zeros(nParamsForHist,1);
for i = 1:nParamsForHist
    qPlusPanels(i) = 3*(i-1) + 1; 
    randPanels(i) = 3*(i-1) + 2;
end

%% Main plot loop through each level of factorName
% Loop through each level of the chosen factorName
for i = 1:length(ind)
    % Grab the next unique level
    level = data.(factorName)(ind(i));
    % Average the parameters for that level
    [qpAverageParams(i,:),randomAverageParams(i,:),qpAvBOLD,randAvBOLD] = oneNoiseLevel(avgResults,factorName,level);
    
    % Select all rows belonging to that level. 
    qpRows = (data.qpPres=="qpControl" & data.(factorName)==level);
    qpRows = data(qpRows==1,:);
    randRows = (data.qpPres=="random" & data.(factorName)==level);
    randRows = data(randRows==1,:);

    % Get the average for qp and random for that level
    averageQPResponse = makePredicted(model, stimDomain, qpAverageParams(i,:), min(stimDomain),qpAvBOLD);
    averageRandomResponse = makePredicted(model, stimDomain, randomAverageParams(i,:), min(stimDomain),randAvBOLD);
    
    % Print the level and the values for the sigma/bold/nOutcomes constant
    % for each level.
    %fprintf('------------------------------------------------------------\n');
    %fprintf('%s: %s\n',factorName,num2str(level));
    %fprintf('maxBOLD / sigma\n');
    %fprintf('Q+ sigma-bold-nOutcomes value =     %.03f\n', mean(qpRows.sbnConstant));
    %fprintf('Random sigma-bold-nOutcomes value = %.03f\n', mean(randRows.sbnConstant));

    
    % Create a new figure for each level
    mainFig = figure;
    set(gcf,'Position',[50 50 1200 700]);
    
    subplot(nParamsForHist,3,qPlusPanels);
    
    hold on;
    for j = 1:size(qpRows,1)
        params = table2array(qpRows(j,1:length(paramNamesInOrder)));
        qpResponse = makePredicted(model, stimDomain, params, min(stimDomain),qpRows.maxBOLD(j));
        plot1 = plot(stimDomain,qpResponse,'-','Color' ,darkQPColor,'LineWidth',2,'HandleVisibility','off');
        plot1.Color(4) = 0.1;
    end
    
    plotV = plot(stimDomain,predictedRelativeResponse,'-','Color',veridicalColor,'LineWidth',6);
    plotV.Color(4) = .75;
    plot(stimDomain,averageQPResponse,'--','Color',darkQPColor,'LineWidth',3);
    
    %ylim([0 max([qpRows.maxBOLD; randRows.maxBOLD])]);
    ylim([0 1]);
    set(gca,'XScale', 'lin');
    xlabel('Contrast');
    ylabel('Predicted Response, Normalized 0-1');
    legend('Veridical','Q+','Location','Northwest');
    title(horzcat(func2str(model), ' curves for ',factorName,': ',num2str(level)));
    hold off;
    
    
    %% Plot the random fits in the middle panel
    subplot(nParamsForHist,3,randPanels);
    hold on;
    
    for j = 1:size(randRows,1)
        params = table2array(randRows(j,1:length(paramNamesInOrder)));
        randomResponse = makePredicted(model, stimDomain, params, min(stimDomain),randRows.maxBOLD(j));
        plot2 = plot(stimDomain,randomResponse,'-','Color',darkRandColor,'LineWidth',2,'HandleVisibility','off');
        plot2.Color(4) = 0.1;
    end
    
    % Plot veridical and average parameter fits
    plotV = plot(stimDomain,predictedRelativeResponse,'-','Color',veridicalColor,'LineWidth',6);
    plotV.Color(4) = .75;
    plot(stimDomain,averageRandomResponse,'--','Color',lightRandColor,'LineWidth',3);
    
    %Formatting
%    ylim([0 max([qpRows.maxBOLD; randRows.maxBOLD])]);
    ylim([0 1]);
    set(gca,'XScale', 'lin');
    xlabel('Contrast');
    ylabel('Predicted Response, Normalized 0-1');
    legend('Veridical','Random','Location','Northwest');
    title(horzcat(func2str(model), ' curves for ',factorName,': ',num2str(level)));
    hold off;
    
    %% Third panel histograms
    
    for j = 1:nParamsForHist
        panel = j*3;
        binRange = paramDomains{j};
        sampleSimulatedToPlot = qpRows.(strcat(paramNamesToPlot{j},'Sim'))(1);
        histogramPanel(paramNamesToPlot{j},nParamsForHist,panel,...,
            binRange,qpRows,randRows,darkQPColor,darkRandColor,veridicalColor,...,
            sampleSimulatedToPlot);
        
    end
   
    % Save main figure.
    set(mainFig,'Units','Inches');
    pos = get(mainFig,'Position');
    set(mainFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(mainFig,horzcat(dirName,'/',func2str(model),'_',num2str(level),'.pdf'),'-dpdf','-r0')

    
    

end

%{

%% Assemble signed and unsigned errors
data.slopeSigned = data.slope - data.slopeSim;
data.semiSatSigned = data.semiSat - data.semiSatSim;
data.betaSigned = data.beta - data.betaSim;
data.sigmaSigned = data.sigma - data.sigmaSim;
data.maxBOLDSigned = data.maxBOLD - data.maxBOLDSim;

data.slopeUnsigned = abs(data.slope - data.slopeSim);
data.semiSatUnsigned = abs(data.semiSat - data.semiSatSim);
data.betaUnsigned = abs(data.beta - data.betaSim);
data.sigmaUnsigned = abs(data.sigma - data.sigmaSim);
data.maxBOLDUnsigned = abs(data.maxBOLD - data.maxBOLDSim);


SignedResults = varfun(@mean,data,'InputVariables',{'slopeSigned','semiSatSigned','betaSigned','sigmaSigned','maxBOLDSigned'},...
    'GroupingVariables',{'qpPres','sigmaSim'})
 
SignedResultsStd = varfun(@std,data,'InputVariables',{'slopeSigned','semiSatSigned','betaSigned','sigmaSigned','maxBOLDSigned'},...
     'GroupingVariables',{'qpPres','sigmaSim'});

UnsignedResults = varfun(@mean,data,'InputVariables',{'slopeUnsigned','semiSatUnsigned','betaUnsigned','sigmaUnsigned','maxBOLDUnsigned'},...
    'GroupingVariables',{'qpPres','sigmaSim'})

UnsignedResultsStd = varfun(@std,data,'InputVariables',{'slopeUnsigned','semiSatUnsigned','betaUnsigned','sigmaUnsigned','maxBOLDUnsigned'},...
     'GroupingVariables',{'qpPres','sigmaSim'});

for i = 1:size(SignedResults,1)
    SignedResults.Level(i) = strcat(SignedResults.qpPres(i),' Noise: ',num2str(SignedResults.sigmaSim(i)));
    UnsignedResults.Level(i) = strcat(UnsignedResults.qpPres(i),' Noise: ',num2str(UnsignedResults.sigmaSim(i)));
end

%Calculate SEM
semSigned = zeros(6,3);
semSigned(:,1) = SignedResultsStd.std_slopeSigned ./ sqrt(SignedResults.GroupCount);
semSigned(:,2) = SignedResultsStd.std_semiSatSigned ./ sqrt(SignedResults.GroupCount);
semSigned(:,3) = SignedResultsStd.std_maxBOLDSigned ./ sqrt(SignedResults.GroupCount);
%Calculate SEM
semUnsigned = zeros(6,3);
semUnsigned(:,1) = UnsignedResultsStd.std_slopeUnsigned ./ sqrt(UnsignedResultsStd.GroupCount);
semUnsigned(:,2) = UnsignedResultsStd.std_semiSatUnsigned ./ sqrt(UnsignedResultsStd.GroupCount);
semUnsigned(:,3) = UnsignedResultsStd.std_maxBOLDUnsigned ./ sqrt(UnsignedResultsStd.GroupCount);


% Grab the data for easier plotting
signedToPlot = SignedResults{1:6, [4:5 8]};
unsignedToPlot = UnsignedResults{1:6, [4:5 8]};

signedFig = figure; 
set(gcf,'Position',[50 50 1200 700]);
bar(signedToPlot, 'grouped');
hold on;
% Find the number of groups and the number of bars in each group
ngroups = size(signedToPlot, 1);
nbars = size(signedToPlot, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, signedToPlot(:,i), semSigned(:,i), 'k', 'linestyle', 'none');
end

legend({'Slope','Semi-Sat','maxBOLD'});
xticklabels(SignedResults{:,9});
ylabel('Signed error (SEM)');
title('Signed Error');

set(signedFig,'Units','Inches');
pos = get(signedFig,'Position');
set(signedFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(signedFig,'./SEM_Signed_Error.pdf','-dpdf','-r0')



unsignedFig = figure; 
set(gcf,'Position',[50 50 1200 700]);
bar(unsignedToPlot, 'grouped');
hold on;
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, unsignedToPlot(:,i), semUnsigned(:,i), 'k', 'linestyle', 'none');
end

legend({'Slope','Semi-Sat','maxBOLD'});
xticklabels(UnsignedResults{:,9});
ylabel('Unsigned error (SEM)');
title('Unsigned Error');
set(unsignedFig,'Units','Inches');
pos = get(unsignedFig,'Position');
set(unsignedFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(unsignedFig,'./SEM_Unsigned_Error.pdf','-dpdf','-r0')



yVeridical = logistic(stimDomain,[sampleSimulatedParams(1) sampleSimulatedParams(2) sampleSimulatedParams(3)]);

yQPLowNoise = logistic(stimDomain,[avgResults.mean_slope(1) avgResults.mean_semiSat(1) 1]);
yRandomLowNoise = logistic(stimDomain,[avgResults.mean_slope(2) avgResults.mean_semiSat(2) 1]);

yQPMedNoise = logistic(stimDomain,[avgResults.mean_slope(3) avgResults.mean_semiSat(3) 1]);
yRandomMedNoise = logistic(stimDomain,[avgResults.mean_slope(4) avgResults.mean_semiSat(4) 1]);

yQPHighNoise = logistic(stimDomain,[avgResults.mean_slope(5) avgResults.mean_semiSat(5) 1]);
yRandomHighNoise = logistic(stimDomain,[avgResults.mean_slope(6) avgResults.mean_semiSat(6) 1]);

fprintf('\n Low Noise  QP: %.04f | Random %.04f',corr(yQPLowNoise',yVeridical'),corr(yRandomLowNoise',yVeridical'));
fprintf('\n Med Noise  QP: %.04f | Random %.04f',corr(yQPMedNoise',yVeridical'),corr(yRandomMedNoise',yVeridical'));
fprintf('\n High Noise QP: %.04f | Random %.04f\n',corr(yQPHighNoise',yVeridical'),corr(yRandomHighNoise',yVeridical'));
%}


%% Useful sub-functions
   
function [qpParams,randomParams,qpAvBOLD,randAvBOLD] = oneNoiseLevel(avgResults,factorName,level)
% Extract average qpParams and randomParams for one factor level
qpRows = (avgResults.qpPres=="qpControl" & avgResults.(factorName)==level);
a = avgResults(qpRows,:);
qpParams = [a.mean_slope a.mean_semiSat a.mean_beta a.mean_sigma];
qpAvBOLD = a.mean_maxBOLD;

randRows = (avgResults.qpPres=="random" & avgResults.(factorName)==level);
b = avgResults(randRows,:);
randomParams = [b.mean_slope b.mean_semiSat b.mean_beta b.mean_sigma];
randAvBOLD = b.mean_maxBOLD;

end

function predictedResponse = makePredicted(model, stimDomain, params, baseline, maxBOLD)
predictedResponse = model(stimDomain,params) - model(baseline,params);
%predictedResponse = predictedResponse .* maxBOLD;
end

function histogramPanel(paramName,nParamsForHist,panel,binRange,qpRows,randRows,darkQPColor,darkRandColor,veridicalColor,sampleSimulatedParam)
    
    hcx = histcounts(qpRows.(paramName),[binRange Inf]);
    hcy = histcounts(randRows.(paramName),[binRange Inf]);
    subplot(nParamsForHist,3,panel);
    plotV = plot([sampleSimulatedParam sampleSimulatedParam],[0 1.25*max([hcx hcy])],'Color',veridicalColor,'LineWidth',4);
    plotV.Color(4) = .75;
    hold on;
    b = bar(binRange,[hcx;hcy]','BarWidth',3);
    b(1).FaceColor = darkQPColor;
    b(2).FaceColor = darkRandColor;
    b(1).FaceAlpha = .8;
    b(2).FaceAlpha = .8;
    ylim([0 1.25*max([hcx hcy])]);
    if panel == 3
        legend('Veridical','Q+','Random');
    end
    xlabel(sprintf('%s Value',paramName));
    ylabel('Number of Simulations');
    title(sprintf('%s parameter estimates',paramName));

    ax = gca;
    xLoc = ax.Position(1) + .22;
    yLoc = ax.Position(2)+.05;
    
    vParam = sprintf('Veridical: %.03f\n',sampleSimulatedParam);
    annotation('textbox',[xLoc yLoc .1 .1],'String',vParam,'EdgeColor','none');
    qPlusParam = sprintf('Q+ M(SD)\n%.03f(%.02f)',mean(qpRows.(paramName)),std(qpRows.(paramName)));
    annotation('textbox',[xLoc yLoc-.04 .1 .1],'String',qPlusParam,'EdgeColor','none');
    randParam = sprintf('Random M(SD)\n%.03f(%.02f)',mean(randRows.(paramName)),std(randRows.(paramName)));
    annotation('textbox',[xLoc yLoc-.1 .1 .1],'String',randParam,'EdgeColor','none');
    
    hold off;
end

end