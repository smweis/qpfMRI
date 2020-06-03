% Plot and print some summary statistics 

%% Provide some info about the simulations
model = @logistic;
sameParams = true; % Are the veridical model params the same for all parameters?
stimDomain = linspace(.01,1,30); % What's the stimDomain?
baseline = .01;
reloadData = true; % Do we need to reload the data into memory?

% What factor do you want to look at? 
%factorName = 'sigmaSim';
factorName = 'nOutcomes';

% How many parameters to look at in the main figure?
paramNamesToPlot = {'slope','semiSat','maxBOLD','sigma'};
paramDomains = {linspace(0,1,100),linspace(0,1,100),linspace(0,2.5,100),linspace(0,8,100)};
nParamsForHist = length(paramNamesToPlot);

%% Check model, grab parameter names, load and sanity check data
paramNamesInOrder = checkModel(model);

% If 'data' not in memory, run the function to extract the data from the directory

if reloadData
    % Parameter set 1: .23 .42; 200 sims
    %[data] = loadSimulatedData('logisticResultsParamSet1/LogisticThreeNoiseLevels',model,[1,2]);
    % Parameter set 2: .42, .79; 100 sims
    %[data] = loadSimulatedData('logisticResultsParamSet2',model,[1,2]);
    % Parameter set 3: .23 .42; 100 sims, changed maxBOLD fit
    %[data] = loadSimulatedData('logisticResultsParamSet3/Results',model,[1,2]);
    %[data] = loadSimulatedData('logisticResultsParamSet4/Results',model,[1,2]);
    %[data] = loadSimulatedData('logisticResultsParamSet5/Results',model,[1,2]);    
    % Parameter set 6: .19, .49, 200 sims
    %[data] = loadSimulatedData('logisticResultsParamSet6/Results',model,[1,2]); 
    % Parameter set 7: .19, .49, sigmaSim = .2; nOutcomes vary: 25, 51, 99
    [data] = loadSimulatedData('logisticResultsParamSet7_nOutcomesTest/Results',model,[1,2]); 
end

%% Check the loaded data for duplicate values. 
% We want to check for duplicate values in our fit parameters. These should
% almost NEVER result in identical parameters; otherwise, we suspect
% something went wrong with setting the random seed. 
[~, dupeInd] = unique(data(:,1:length(paramNamesInOrder)+1), 'rows');
duplicate_ind = setdiff(1:size(data, 1:length(paramNamesInOrder)+1), dupeInd);
if ~isempty(duplicate_ind)
    warning('Duplicate values detected. Inspect data and logs to ensure unique simulations.');
end

%% Print averages by Q+ selection and a factor of your choice

% If params are all the same, we can use this: 
if sameParams
    sampleSimulatedParams = struct;
    sampleSimulatedParamsVec = zeros(1,length(paramNamesInOrder));
    fprintf('One veridical parameter set is used:\n');
    for j = 1:nParamsForHist
        sampleSimulatedParams.(paramNamesToPlot{j}) = data.(strcat(paramNamesToPlot{j},'Sim'))(1);
        fprintf('Simulated parameter "%s" = %.03\n',paramNamesToPlot{j},data.(strcat(paramNamesToPlot{j},'Sim'))(1));
    end
    for j = 1:length(paramNamesInOrder)
        sampleSimulatedParamsVec(j) = data.(strcat(paramNamesInOrder{j},'Sim'))(1);
    end
    predictedRelativeResponse = makePredicted(model, stimDomain, sampleSimulatedParamsVec, baseline);
else
    error('multiple params not handled at this point');
end

% Print the average result for Q+/random and the given factor. 
avgResults = varfun(@mean,data,'InputVariables',[paramNamesInOrder 'maxBOLD'],...
       'GroupingVariables',{factorName,'qpPres'})



%% Plotting colors
% Enter Colors for plot
lightRandColor = '#9E91FF';
darkRandColor = '#9E00FF';
lightQPColor = '#4CB963';
darkQPColor = '#008148';
veridicalColor = '#FC7A1E';



% Detect number of levels for the given factor name:
[~, ind] = unique(data.(factorName), 'rows');
qpAverageParams = zeros(length(ind),length(paramNamesInOrder));
randomAverageParams = zeros(length(ind),length(paramNamesInOrder));

% Finding the relation between sigma, maxBOLD, and nOutcomes
data.sbnConstant = (data.maxBOLD ./ data.sigma);


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
    [qpAverageParams(i,:),randomAverageParams(i,:)] = oneNoiseLevel(avgResults,factorName,level);
    
    % Select all rows belonging to that level. 
    qpRows = (data.qpPres=="qpControl" & data.(factorName)==level);
    qpRows = data(qpRows==1,:);
    randRows = (data.qpPres=="random" & data.(factorName)==level);
    randRows = data(randRows==1,:);

    % Get the average for qp and random for that level
    averageQPResponse = makePredicted(model, stimDomain, qpAverageParams(i,:), min(stimDomain));
    averageRandomResponse = makePredicted(model, stimDomain, randomAverageParams(i,:), min(stimDomain));
    
    % Print the level and the values for the sigma/bold/nOutcomes constant
    % for each level.
    fprintf('------------------------------------------------------------\n');
    fprintf('%s: %s\n',factorName,num2str(level));
    fprintf('maxBOLD / sigma\n');
    fprintf('Q+ sigma-bold-nOutcomes value =     %.03f\n', mean(qpRows.sbnConstant));
    fprintf('Random sigma-bold-nOutcomes value = %.03f\n', mean(randRows.sbnConstant));

    
    % Create a new figure for each level
    mainFig = figure;
    set(gcf,'Position',[50 50 1200 700]);
    
    subplot(nParamsForHist,3,qPlusPanels);
    
    hold on;
    for j = 1:size(qpRows,1)
        params = table2array(qpRows(j,1:length(paramNamesInOrder)));
        qpResponse = makePredicted(model, stimDomain, params, min(stimDomain));
        plot1 = plot(stimDomain,qpResponse,'-','Color' ,lightQPColor,'LineWidth',.8,'HandleVisibility','off');
        plot1.Color(4) = 0.1;
    end
    

    plot(stimDomain,averageQPResponse,':','Color',darkQPColor,'LineWidth',4);
    plotV = plot(stimDomain,predictedRelativeResponse,'-','Color',veridicalColor,'LineWidth',6);
    plotV.Color(4) = .5;
    ylim([0 1]);
    set(gca,'XScale', 'lin');
    xlabel('Contrast');
    ylabel('Predicted Response, Normalized 0-1');
    legend('Q+','Veridical','Location','Northwest');
    title(horzcat(func2str(model), ' Curves for ',factorName,': ',num2str(level)));
    hold off;
    
    
    %% Plot the random fits in the middle panel
    subplot(nParamsForHist,3,randPanels);
    hold on;
    for j = 1:size(randRows,1)
        params = table2array(randRows(j,1:length(paramNamesInOrder)));
        randomResponse = makePredicted(model, stimDomain, params, min(stimDomain));
        plot2 = plot(stimDomain,randomResponse,'-','Color',lightRandColor,'LineWidth',.8,'HandleVisibility','off');
        plot2.Color(4) = 0.1;
    end
    % Plot veridical and average parameter fits
    plot(stimDomain,averageRandomResponse,':','Color',darkRandColor,'LineWidth',4);
    plotV = plot(stimDomain,predictedRelativeResponse,'-','Color',veridicalColor,'LineWidth',6);
    plotV.Color(4) = .5;
    
    %Formatting
    ylim([0 1]);
    set(gca,'XScale', 'lin');
    xlabel('Contrast');
    ylabel('Predicted Response, Normalized 0-1');
    legend('Random','Veridical','Location','Northwest');
    title(horzcat(func2str(model), ' Curves for ',factorName,': ',num2str(level)));
    hold off;
    
    %% Third panel histograms
    
    for j = 1:nParamsForHist
        panel = j*3;
        binRange = paramDomains{j};
        histogramPanel(paramNamesToPlot{j},nParamsForHist,panel,...,
            binRange,qpRows,randRows,darkQPColor,darkRandColor,veridicalColor,...,
            sampleSimulatedParams.(paramNamesToPlot{j}));
        
    end
   
    % Save main figure.
    set(mainFig,'Units','Inches');
    pos = get(mainFig,'Position');
    set(mainFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(mainFig,horzcat('./',func2str(model),'_',num2str(level),'.pdf'),'-dpdf','-r0')

    
    

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
   
function [qpParams,randomParams] = oneNoiseLevel(avgResults,factorName,level)
% Extract average qpParams and randomParams for one factor level
qpRows = (avgResults.qpPres=="qpControl" & avgResults.(factorName)==level);
a = avgResults(qpRows,:);
qpParams = [a.mean_slope a.mean_semiSat a.mean_beta a.mean_sigma];

randRows = (avgResults.qpPres=="random" & avgResults.(factorName)==level);
b = avgResults(randRows,:);
randomParams = [b.mean_slope b.mean_semiSat b.mean_beta b.mean_sigma];


end

function predictedResponse = makePredicted(model, stimDomain, params, baseline)
predictedResponse = model(stimDomain,params) - model(baseline,params);
end

function histogramPanel(paramName,nParamsForHist,panel,binRange,qpRows,randRows,darkQPColor,darkRandColor,veridicalColor,sampleSimulatedParam)
    
    hcx = histcounts(qpRows.(paramName),[binRange Inf]);
    hcy = histcounts(randRows.(paramName),[binRange Inf]);
    subplot(nParamsForHist,3,panel);
    b = bar(binRange,[hcx;hcy]','BarWidth',1.5);
    b(1).FaceColor = darkQPColor;
    b(2).FaceColor = darkRandColor;
    hold on;
    plot([sampleSimulatedParam sampleSimulatedParam],ylim,'Color',veridicalColor,'LineWidth',2)
    legend('Q+','Random','Veridical');
    xlabel(sprintf('%s Value',paramName));
    ylabel('Number of Simulations');
    title(sprintf('Histogram of %s Parameter Estimates',paramName));
    
    numPlot = panel/3;
    yLoc = .8 - (numPlot-1)*(.6/(nParamsForHist - 1));
    
    vParam = sprintf('Veridical: %.03f\n',sampleSimulatedParam);
    annotation('textbox',[.9 yLoc .1 .1],'String',vParam,'EdgeColor','none');
    qPlusParam = sprintf('Q+ M(SD)\n%.03f(%.02f)',mean(qpRows.(paramName)),std(qpRows.(paramName)));
    annotation('textbox',[.9 yLoc-.02 .1 .1],'String',qPlusParam,'EdgeColor','none');
    randParam = sprintf('Random M(SD)\n%.03f(%.02f)',mean(randRows.(paramName)),std(randRows.(paramName)));
    annotation('textbox',[.9 yLoc-.1 .1 .1],'String',randParam,'EdgeColor','none');
    
    hold off;
end