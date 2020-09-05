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
paramsDomain = struct;
paramsDomain.slope = makeDomain(-1.5,-.2,10,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.maxBOLD = makeDomain(.75,1.25,11,'spacing','zeno');
paramsDomain.sigma = makeDomain(.5,4,8);
nParamsForHist = 0;
paramNamesToPlot = {};
factorName = 'nOutcomes';
sameParams = true; 
stimDomain = makeDomain(.01,1,25);
baseline = .01;
posterMode = true;
dirStem = pwd;
dirName = fullfile(dirStem,'results','7vs31_.15sdNoise');

data = summarizeAndPlotSimulations(model,paramsDomain,factorName,...,
    sameParams,stimDomain,baseline,dirName,'posterMode',posterMode,...,
    'nParamsForHist',nParamsForHist,'paramNamesToPlot',paramNamesToPlot);


% Run again without re-loading data
data = summarizeAndPlotSimulations(model,paramsDomain,factorName,...,
    sameParams,stimDomain,baseline,dirName,data,'posterMode',posterMode,...,
    'nParamsForHist',nParamsForHist,'paramNamesToPlot',paramNamesToPlot);

---------------------------------------------------------------------------
% Example 2: Make figures for all parameter sets: 
model = @logistic;
paramsDomain = struct;
paramsDomain.slope = makeDomain(-1.2,-.2,20,'spacing','log');
paramsDomain.semiSat = makeDomain(.01,1,10);
paramsDomain.maxBOLD = makeDomain(0,2.5,100);
paramsDomain.sigma = makeDomain(0,9,100);
nParamsForHist = 0;
paramNamesToPlot = {};
factorName = 'sigmaSim';
sameParams = true; 
stimDomain = makeDomain(.01,1,25);
baseline = .01;
posterMode = true;

a = dir('logistic*');

for name = 1:length(a) 
    data = summarizeAndPlotSimulations(model,paramsDomain,factorName,...,
    sameParams,stimDomain,baseline,a(name).name,'posterMode',posterMode,...,
    'nParamsForHist',nParamsForHist,'paramNamesToPlot',paramNamesToPlot);
    close all; 
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
p.addParameter('nParamsForHist',0,@isnumeric);
p.addParameter('paramNamesToPlot',{},@iscell);

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
    [data] = loadSimulatedData(dirName,model,[1,2]);
else
    fprintf('Data being re-analyzed.');
    data = p.Results.data;
end


%% Print averages by Q+ selection and a factor of your choice
assert(length(p.Results.paramNamesToPlot) == p.Results.nParamsForHist,'paramsDomain should only include the parameters you want to see plotted as histograms');
    
% If params are all the same, we can use this: 
if sameParams
    sampleSimulatedParams = struct;
    sampleSimulatedParamsVec = zeros(1,length(paramNamesInOrder));
    fprintf('\nOne veridical parameter set is used:\n');
    for j = 1:p.Results.nParamsForHist
        sampleSimulatedParams.(p.Results.paramNamesToPlot{j}) = data.(strcat(p.Results.paramNamesToPlot{j},'Sim'))(1);
        fprintf('Simulated parameter "%s" = %.03f\n',p.Results.paramNamesToPlot{j},data.(strcat(p.Results.paramNamesToPlot{j},'Sim'))(1));
    end
    for j = 1:length(paramNamesInOrder)
        sampleSimulatedParamsVec(j) = data.(strcat(paramNamesInOrder{j},'Sim'))(1);
        sampleMaxBOLD = data.maxBOLDSim(1);
    end
    modelResponseNoCorrection = model(stimDomain,sampleSimulatedParamsVec);
    predictedRelativeResponse = makePredicted(model, stimDomain, sampleSimulatedParamsVec, baseline, sampleMaxBOLD);
else
    error('multiple params not handled at this point');
end


% Print the average result for Q+/random and the given factor. 
avgResults = varfun(p.Results.avgFunc,data,'InputVariables',[paramNamesInOrder 'maxBOLD' 'nOutcomes'],...
       'GroupingVariables',{factorName,'qpPres'})

avgFuncName = func2str(p.Results.avgFunc);

%% Plotting colors



% Detect number of levels for the given factor name:
[~, ind] = unique(data.(factorName), 'rows');
qpAverageParams = zeros(length(ind),length(paramNamesInOrder));
randomAverageParams = zeros(length(ind),length(paramNamesInOrder));

% Finding the relation between sigma, maxBOLD, and nOutcomes
%data.sbnConstant = (data.maxBOLD ./ data.sigma);

if p.Results.nParamsForHist > 0
    subPlotRows = p.Results.nParamsForHist;
    subPlotCols = 3;
    qPlusPanels = zeros(p.Results.nParamsForHist,1);
    randPanels = zeros(p.Results.nParamsForHist,1);
    for i = 1:p.Results.nParamsForHist
        qPlusPanels(i) = 3*(i-1) + 1; 
        randPanels(i) = 3*(i-1) + 2;
    end
else
    subPlotRows = 2;
    subPlotCols = 2;
    qPlusPanels = [1 3];
    randPanels = [2 4];
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
    set(gcf,'color','w');
    subplot(subPlotRows,subPlotCols,qPlusPanels);
    
        
    hold on;
    for j = 1:size(qpRows,1)
        params = table2array(qpRows(j,1:length(paramNamesInOrder)));
        qpResponse = makePredicted(model, stimDomain, params, min(stimDomain),qpRows.maxBOLD(j));
        simResponse = model(stimDomain,params);
        qpRows.rSquared(j) = corr(modelResponseNoCorrection',simResponse')^2;
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
    else
        plotFactorName = factorName;
    end
    
    sgtitle(sprintf('%s Model Fits for %s = %s(SD)',modelName,plotFactorName,num2str(level)),'FontSize',45);
    hold off;
    
    
    %% Plot the random fits in the middle panel
    subplot(subPlotRows,subPlotCols,randPanels);
    hold on;
    
    for j = 1:size(randRows,1)
        params = table2array(randRows(j,1:length(paramNamesInOrder)));
        randomResponse = makePredicted(model, stimDomain, params, min(stimDomain),randRows.maxBOLD(j));
        simResponse = model(stimDomain,params);
        randRows.rSquared(j) = corr(modelResponseNoCorrection',simResponse')^2;
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
    
    for j = 1:p.Results.nParamsForHist
        panel = j*3;
        binRange = paramsDomain.(p.Results.paramNamesToPlot{j});
        sampleSimulatedToPlot = qpRows.(strcat(p.Results.paramNamesToPlot{j},'Sim'))(1);
        histogramPanel(p.Results.paramNamesToPlot{j},p.Results.nParamsForHist,panel,...,
            binRange,qpRows,randRows,sampleSimulatedToPlot,p.Results.avgFunc,posterFormat);
        
    end
    paramString = '';
    for m = 1:length(paramNamesInOrder)-2
        paramString = [paramString newline sprintf('%s = %.02f',paramNamesInOrder{m},sampleSimulatedParamsVec(m))];
    end
    simParamString = sprintf('Simulated Parameters:\n%s',paramString);
    annotation('textbox',[.2 .7 .1 .1],'String',simParamString,'EdgeColor','none','FontSize',15);
    
    % Save main figure.

    set(mainFig,'Units','Inches');
    pos = get(mainFig,'Position');
    set(mainFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(mainFig,horzcat(dirName,'/',func2str(model),'_',num2str(level),'.pdf'),'-dpdf','-r0')

    
    screen = get(0, 'Screensize');
    screen(3) = screen(3) - 400;
    rFig = figure('Position',screen);
    set(gcf,'color','w');
    hold on;
    hcx = histcounts(qpRows.rSquared,[linspace(.8,1,40) Inf]);
    hcy = histcounts(randRows.rSquared,[linspace(.8,1,40) Inf]);
    b = bar(linspace(.8,1,40),[hcx;hcy]','BarWidth',3);
    b(1).FaceColor = posterFormat.darkQPColor;
    b(2).FaceColor = posterFormat.darkRandColor;
    b(1).FaceAlpha = .8;
    b(2).FaceAlpha = .8;
    [h,pVal,~,stats] = ttest2(qpRows.rSquared,randRows.rSquared);
    ax = gca;
    ax.FontSize = 40;
    xlabel('R^2 Values','FontSize',60);
    ylabel('Number of Simulations','FontSize',60);
    %title(sprintf('R^2 %s Level = %s',plotFactorName,num2str(level)),'FontSize',45);
    
    xLoc = ax.Position(1) + .1;
    yLoc = ax.Position(2) + .4;
    
    if h == 0
        ttestString = sprintf('No Significant Difference\n t(%d) = %.02f, p = %.05f\nQ+ Median = %.03f(%.03f)\nRand Median = %.03f(%.03f)',...,
            stats.df,stats.tstat,pVal,median(qpRows.rSquared),std(qpRows.rSquared),median(randRows.rSquared),std(randRows.rSquared));
    else
        ttestString = sprintf('t(%d) = %.02f, p = %.05f\nQ+ Median = %.03f(%.03f)\nRand Median = %.03f(%.03f)',...,
            stats.df,stats.tstat,pVal,median(qpRows.rSquared),std(qpRows.rSquared),median(randRows.rSquared),std(randRows.rSquared));
    end
    annotation('textbox',[xLoc yLoc .1 .1],'String',ttestString,'EdgeColor','none','FontSize',15);
    

    annotation('textbox',[xLoc yLoc - .2 .1 .1],'String',simParamString,'EdgeColor','none','FontSize',15);
    
    set(rFig,'Units','Inches');
    pos = get(rFig,'Position');
    set(rFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(rFig,horzcat(dirName,'/',func2str(model),'_',num2str(level),'rSquared.pdf'),'-dpdf','-r0')

    hold off;

end


permTests = 1000;
permResults = zeros(1,permTests);
randParams = zeros(1,length(paramNamesInOrder)-1);
for i = 1:permTests
    for j = 1:length(paramNamesInOrder)-2
        randParams(j) = randsample(paramsDomain.(paramNamesInOrder{j}),1);
    end
    randParams(end) = 1;
    randResponse = model(stimDomain,randParams);
    permResults(i) = corr(randResponse',modelResponseNoCorrection')^2;
end





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
%predictedResponse = model(stimDomain,params);
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
    
    if nParamsForHist == 2
        yLoc1 = ax.Position(2) + .1;
        yLoc2 = yLoc1 - .07;
    elseif nParamsForHist == 3
        yLoc1 = ax.Position(2) + .05;
        yLoc2 = yLoc1 - .07;
    else
        yLoc1 = ax.Position(2);
        yLoc2 = yLoc1 - .07;
    end
    
    qPlusParam = sprintf('%s(SD)\n%.02f(%.02f)',func2str(avgFunc),avgFunc(qpRows.(paramName)),std(qpRows.(paramName)));
    annotation('textbox',[xLoc yLoc1 .1 .1],'String',qPlusParam,'EdgeColor','none','FontSize',15,'Color',posterFormat.darkQPColor);
    randParam = sprintf('%s(SD)\n%.02f(%.02f)',func2str(avgFunc),avgFunc(randRows.(paramName)),std(randRows.(paramName)));
    annotation('textbox',[xLoc yLoc2 .1 .1],'String',randParam,'EdgeColor','none','FontSize',15,'Color',posterFormat.darkRandColor);
    
    hold off;
end

end