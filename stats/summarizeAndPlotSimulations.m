function [data]= summarizeAndPlotSimulations(model,paramsDomain,factorName,sameParams,stimDomain,baseline,dirName,varargin)
%% [data]= summarizeAndPlotSimulations(model,paramsDomain,factorName,sameParams,stimDomain,baseline,dirName,varargin)

% Plot and print some summary statistics for simulations results


% Syntax:
%  [data]= summarizeAndPlotSimulations(model,paramsDomain,factorName,sameParams,stimDomain,baseline,dirName,varargin)
%
% Description:
%	Plots summary statistics for a set of simulations.
%
% Inputs:
%   model                 - A function handle. This should be the
%                           'continuous function.' The Quest+ specific
%                           function will be defined in the model-specific
%                           code block. Currently supported:
%                             @doeTemporalModel
%                             @watsonTemporalModel
%                             @logistic
%   paramsDomain          - Struct. Domains for all parameters to plot.
%   factorName            - String. Name of value that varied (will create
%                           separate plots for each value.)
%   sameParams            - Boolean. Were all simulations using the same
%                           parameter values (with the exception of factorName)?
%   stimDomain            - Vector. Values that the stimulus were drawn
%                           from.
%   baseline              - Numeric. Value to be used for the baseline
%                           trial.
%   dirName               - String. Full name of the directory to find the
%                           data. 
% Optional positional input
%   data                  - Table. The loaded-in data. This is the high
%                           time-cost of this function, so if you'd like to 
%                           re-load the data, add it as the last optional argument.
% Optional key/value pairs   
%   'avgFunc'              - Function handle. Default (@median). Function
%                            to use to evaluate the average parameters.
%   'posterMode'           - Boolean. Default = false. Will alter some
%                            formatting if we're plotting for posters. 
%{
---------------------------------------------------------------------------
% Example 1:

model = @logistic;

% Here, only include parameters to plot.
paramsDomain = struct;
paramsDomain.slope = makeDomain(.01,1,100);
paramsDomain.semiSat = makeDomain(.01,1,100);
%paramsDomain.maxBOLD = makeDomain(0,2.5,100);
%paramsDomain.sigma = makeDomain(0,9,100);
factorName = 'sigmaSim';
sameParams = true; 
stimDomain = makeDomain(.01,1,25);
baseline = .01;
posterMode = true;
dirStem = pwd;
dirName = fullfile(dirStem,'logisticResultsParamSet9');

data = summarizeAndPlotSimulations(model,paramsDomain,factorName,...,
    sameParams,stimDomain,baseline,dirName,'posterMode',posterMode);


% Run again without re-loading data
data = summarizeAndPlotSimulations(model,paramsDomain,factorName,...,
    sameParams,stimDomain,baseline,dirName,data,'posterMode',posterMode);
---------------------------------------------------------------------------
% Example 2: Make figures for all parameter sets: 

model = @logistic;

% Here, only include parameters to plot.
paramsDomain = struct;
paramsDomain.slope = makeDomain(.01,1,100);
paramsDomain.semiSat = makeDomain(.01,1,100);
paramsDomain.maxBOLD = makeDomain(0,2.5,100);

factorName = 'sigmaSim';
sameParams = true; 
stimDomain = makeDomain(.01,1,25);
baseline = .01;


a = dir('logistic*');

for name = 1:length(a) 
    data = summarizeAndPlotSimulations(model,paramsDomain,factorName,...,
    sameParams,stimDomain,baseline,a(name).name);
end
%}

%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));
p.addRequired('paramsDomain',@isstruct);
p.addRequired('factorName',@isstr);
p.addRequired('sameParams',@islogical);
p.addRequired('stimDomain',@isvector);
p.addRequired('baseline',@isnumeric);
p.addRequired('dirName',@isstr);

% Optional params
p.addOptional('data',table,@istable);
p.addParameter('avgFunc',@median,@(x) isa(x,'function_handle'));
p.addParameter('posterMode',false,@islogical);
% Parse

p.parse(model,paramsDomain,factorName,sameParams,stimDomain,baseline,dirName,varargin{:});

%% Check model, grab parameter names, load and sanity check da\ta
paramNamesInOrder = checkModel(model);


%% Some formatting for poster mode

posterFormat = struct;
posterFormat.posterMode = p.Results.posterMode;
posterFormat.lightRandColor = '#3D92C9';
posterFormat.darkRandColor = '#165172';
posterFormat.lightQPColor = '#F99C16';
posterFormat.darkQPColor = '#FA4515';
posterFormat.veridicalColor = '#000000';


% If 'data' not in memory, run the function to extract the data from the directory
if isempty(p.Results.data)
    fprintf('Data being re-loaded from %s',dirName)
    fullDirName = fullfile(dirName,'Results');
    [data] = loadSimulatedData(fullDirName,model,[1,2]);
else
    fprintf('Data being re-analyzed.');
    data = p.Results.data;
end


%% Print averages by Q+ selection and a factor of your choice
paramNamesToPlot = fieldnames(paramsDomain);
nParamsForHist = length(paramNamesToPlot);
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
avgResults = varfun(p.Results.avgFunc,data,'InputVariables',[paramNamesInOrder 'maxBOLD'],...
       'GroupingVariables',{factorName,'qpPres'})

avgFuncName = func2str(p.Results.avgFunc);

%% Plotting colors



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
    [qpAverageParams(i,:),randomAverageParams(i,:),qpAvBOLD,randAvBOLD] = oneNoiseLevel(avgResults,factorName,level,avgFuncName);
    
    % Select all rows belonging to that level. 
    qpRows = (data.qpPres=="qpControl" & data.(factorName)==level);
    qpRows = data(qpRows==1,:);
    randRows = (data.qpPres=="random" & data.(factorName)==level);
    randRows = data(randRows==1,:);

    % Get the average for qp and random for that level
    averageQPResponse = makePredicted(model, stimDomain, qpAverageParams(i,:), min(stimDomain),qpAvBOLD);
    averageRandomResponse = makePredicted(model, stimDomain, randomAverageParams(i,:), min(stimDomain),randAvBOLD);

    
    % Create a new figure for each level
    mainFig = figure('Position', get(0, 'Screensize'));
    
    subplot(nParamsForHist,3,qPlusPanels);
    
    hold on;
    for j = 1:size(qpRows,1)
        params = table2array(qpRows(j,1:length(paramNamesInOrder)));
        qpResponse = makePredicted(model, stimDomain, params, min(stimDomain),qpRows.maxBOLD(j));
        plot1 = plot(stimDomain,qpResponse,'-','Color' ,posterFormat.darkQPColor,'LineWidth',2,'HandleVisibility','off');
        plot1.Color(4) = 0.1;
    end
    
    plotV = plot(stimDomain,predictedRelativeResponse,'-','Color',posterFormat.veridicalColor,'LineWidth',6);
    plotV.Color(4) = .75;
    plot(stimDomain,averageQPResponse,'--','Color',posterFormat.darkQPColor,'LineWidth',3);
    
    ylim([0 1]);
    ax = gca;
    ax.FontSize = 25;
    set(gca,'XScale', 'lin');
    xlabel('Contrast','FontSize',30);
    ylabel('Normalized Predicted Response','FontSize',30);
    %legend('Veridical','Q+','Location','Northwest');
    modelName = func2str(model);
    modelName(1) = upper(modelName(1));
    if strcmpi(factorName,'sigmaSim')
        plotFactorName = 'Noise';
    end
    sgtitle(sprintf('%s Model Fits for %s = %s(SD)',modelName,plotFactorName,num2str(level)),'FontSize',45);
    hold off;
    
    
    %% Plot the random fits in the middle panel
    subplot(nParamsForHist,3,randPanels);
    hold on;
    
    for j = 1:size(randRows,1)
        params = table2array(randRows(j,1:length(paramNamesInOrder)));
        randomResponse = makePredicted(model, stimDomain, params, min(stimDomain),randRows.maxBOLD(j));
        plot2 = plot(stimDomain,randomResponse,'-','Color',posterFormat.darkRandColor,'LineWidth',2,'HandleVisibility','off');
        plot2.Color(4) = 0.1;
    end
    
    % Plot veridical and average parameter fits
    plotV = plot(stimDomain,predictedRelativeResponse,'-','Color',posterFormat.veridicalColor,'LineWidth',6);
    plotV.Color(4) = .75;
    plot(stimDomain,averageRandomResponse,'--','Color',posterFormat.lightRandColor,'LineWidth',3);
    
    %Formatting
    ylim([0 1]);
    ax = gca;
    ax.FontSize = 25;
    set(gca,'XScale', 'lin');
    xlabel('Contrast','FontSize',30);
    
    %legend('Veridical','Random','Location','Northwest');
    
    hold off;
    
    %% Third panel histograms
    
    for j = 1:nParamsForHist
        panel = j*3;
        binRange = paramsDomain.(paramNamesToPlot{j});
        sampleSimulatedToPlot = qpRows.(strcat(paramNamesToPlot{j},'Sim'))(1);
        histogramPanel(paramNamesToPlot{j},nParamsForHist,panel,...,
            binRange,qpRows,randRows,sampleSimulatedToPlot,p.Results.avgFunc,posterFormat);
        
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


SignedResults = varfun(p.Results.avgFunc,data,'InputVariables',{'slopeSigned','semiSatSigned','betaSigned','sigmaSigned','maxBOLDSigned'},...
    'GroupingVariables',{'qpPres','sigmaSim'})
 
SignedResultsStd = varfun(@std,data,'InputVariables',{'slopeSigned','semiSatSigned','betaSigned','sigmaSigned','maxBOLDSigned'},...
     'GroupingVariables',{'qpPres','sigmaSim'});

UnsignedResults = varfun(p.Results.avgFunc,data,'InputVariables',{'slopeUnsigned','semiSatUnsigned','betaUnsigned','sigmaUnsigned','maxBOLDUnsigned'},...
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

yQPLowNoise = logistic(stimDomain,[avgResults.(strcat(avgFuncName,'_slope))(1) avgResults.(strcat(avgFuncName,'_semiSat(1) 1]);
yRandomLowNoise = logistic(stimDomain,[avgResults.(strcat(avgFuncName,'_slope(2) avgResults.(strcat(avgFuncName,'_semiSat(2) 1]);

yQPMedNoise = logistic(stimDomain,[avgResults.(strcat(avgFuncName,'_slope(3) avgResults.(strcat(avgFuncName,'_semiSat(3) 1]);
yRandomMedNoise = logistic(stimDomain,[avgResults.(strcat(avgFuncName,'_slope(4) avgResults.(strcat(avgFuncName,'_semiSat(4) 1]);

yQPHighNoise = logistic(stimDomain,[avgResults.(strcat(avgFuncName,'_slope(5) avgResults.(strcat(avgFuncName,'_semiSat(5) 1]);
yRandomHighNoise = logistic(stimDomain,[avgResults.(strcat(avgFuncName,'_slope(6) avgResults.(strcat(avgFuncName,'_semiSat(6) 1]);

fprintf('\n Low Noise  QP: %.04f | Random %.04f',corr(yQPLowNoise',yVeridical'),corr(yRandomLowNoise',yVeridical'));
fprintf('\n Med Noise  QP: %.04f | Random %.04f',corr(yQPMedNoise',yVeridical'),corr(yRandomMedNoise',yVeridical'));
fprintf('\n High Noise QP: %.04f | Random %.04f\n',corr(yQPHighNoise',yVeridical'),corr(yRandomHighNoise',yVeridical'));
%}


%% Useful sub-functions
   
function [qpParams,randomParams,qpAvBOLD,randAvBOLD] = oneNoiseLevel(avgResults,factorName,level,avgFuncName)
% Extract average qpParams and randomParams for one factor level
qpRows = (avgResults.qpPres=="qpControl" & avgResults.(factorName)==level);
a = avgResults(qpRows,:);
qpParams = [a.(strcat(avgFuncName,'_slope')) a.(strcat(avgFuncName,'_semiSat')) a.(strcat(avgFuncName,'_beta')) a.(strcat(avgFuncName,'_sigma'))];
qpAvBOLD = a.(strcat(avgFuncName,'_maxBOLD'));

randRows = (avgResults.qpPres=="random" & avgResults.(factorName)==level);
b = avgResults(randRows,:);
randomParams = [b.(strcat(avgFuncName,'_slope')) b.(strcat(avgFuncName,'_semiSat')) b.(strcat(avgFuncName,'_beta')) b.(strcat(avgFuncName,'_sigma'))];
randAvBOLD = b.(strcat(avgFuncName,'_maxBOLD'));

end

function predictedResponse = makePredicted(model, stimDomain, params, baseline, maxBOLD)
predictedResponse = model(stimDomain,params) - model(baseline,params);
%predictedResponse = predictedResponse .* maxBOLD;
end

function histogramPanel(paramName,nParamsForHist,panel,binRange,qpRows,randRows,sampleSimulatedParam,avgFunc,posterFormat)
    
    hcx = histcounts(qpRows.(paramName),[binRange Inf]);
    hcy = histcounts(randRows.(paramName),[binRange Inf]);
    subplot(nParamsForHist,3,panel);
    plotV = plot([sampleSimulatedParam sampleSimulatedParam],[0 1.25*max([hcx hcy])],'Color',posterFormat.veridicalColor,'LineWidth',4);
    plotV.Color(4) = .4;
    hold on;
    b = bar(binRange,[hcx;hcy]','BarWidth',3);
    b(1).FaceColor = posterFormat.darkQPColor;
    b(2).FaceColor = posterFormat.darkRandColor;
    b(1).FaceAlpha = .8;
    b(2).FaceAlpha = .8;
    ylim([0 1.25*max([hcx hcy])]);
    if panel == 3 && ~posterFormat.posterMode
        legend('Veridical','Q+','Random','FontSize',12,'Orientation','horizontal');
    end
    xlabel('');
    ylabel('# of Simulations','FontSize',12);
    title(sprintf('%s: Veridical = %.02f',paramName,sampleSimulatedParam),'FontSize',25);

    set(gca,'box','off');
    ax = gca;
    ax.FontSize = 15;
    xLoc = ax.Position(1) + .22;
    yLoc = ax.Position(2)+.075;

    qPlusParam = sprintf('%s(SD)\n%.02f(%.02f)',func2str(avgFunc),avgFunc(qpRows.(paramName)),std(qpRows.(paramName)));
    annotation('textbox',[xLoc yLoc .1 .1],'String',qPlusParam,'EdgeColor','none','FontSize',15,'Color',posterFormat.darkQPColor);
    randParam = sprintf('%s(SD)\n%.02f(%.02f)',func2str(avgFunc),avgFunc(randRows.(paramName)),std(randRows.(paramName)));
    annotation('textbox',[xLoc yLoc-.1 .1 .1],'String',randParam,'EdgeColor','none','FontSize',15,'Color',posterFormat.darkRandColor);
    
    hold off;
end

end