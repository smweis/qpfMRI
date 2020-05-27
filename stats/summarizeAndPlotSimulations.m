% Plot and print some summary statistics 

%% Provide some info about the simulations
model = @logistic;
sameParams = true; % Are the veridical params the same for all parameters?
stimDomain = linspace(.01,1,30); % What's the stimDomain?
baseline = .01;
%% Check model, grab parameter names, load and sanity check data
paramNamesInOrder = checkModel(model);

% If 'data' not in memory, run the function to extract the data from the directory

if ~exist('data','var') || size(data,1) == 0
    [data] = loadSimulatedData('logisticResultsParamSet2',model,[1,2]);
end

% We want to check for duplicate values in our fit parameters. These should
% almost NEVER result in identical parameters; otherwise, we suspect
% something went wrong with setting the random seed. 
[~, ind] = unique(data(:,1:length(paramNamesInOrder)+1), 'rows');
duplicate_ind = setdiff(1:size(data, 1:length(paramNamesInOrder)+1), ind);
if ~isempty(duplicate_ind)
    warning('Duplicate values detected. Inspect data and logs to ensure unique simulations.');
end

%% Print averages by Q+ and noise level
avgResults = varfun(@mean,data,'InputVariables',paramNamesInOrder,...
       'GroupingVariables',{'sigmaSim','qpPres'})

  

% If params are all the same, we can use this: 
if sameParams
    sampleSimulatedParams = [data.slopeSim(1) data.semiSatSim(1) data.betaSim(1),data.sigmaSim(1)]
    predictedRelativeResponse = makePredicted(model, stimDomain, sampleSimulatedParams, baseline);
else
    error('multiple params not handled at this point');
end

% Detect number of noise levels:
[~, ind] = unique(data.sigmaSim, 'rows');
qpAverageParams = zeros(length(ind),length(paramNamesInOrder));
randomAverageParams = zeros(length(ind),length(paramNamesInOrder));

% Enter Colors for plot
lightRandColor = '#9E91FF';
darkRandColor = '#9E00FF';
lightQPColor = '#4CB963';
darkQPColor = '#008148';
veridicalColor = '#FC7A1E';

% Loop through each noiseLevel
for i = 1:length(ind)
    % Grab the next unique noise level
    noiseLevel = data.sigmaSim(ind(i));
    % Average the parameters for that level
    [qpAverageParams(i,:),randomAverageParams(i,:)] = oneNoiseLevel(avgResults,noiseLevel);
    
    % Select all rows belonging to that level. 
    qpRows = (data.qpPres=="qpControl" & data.sigmaSim==noiseLevel);
    qpRows = data(qpRows==1,:);
    randRows = (data.qpPres=="random" & data.sigmaSim==noiseLevel);
    randRows = data(randRows==1,:);

    % Get the average for qp and random for that noise level
    averageQPResponse = makePredicted(model, stimDomain, qpAverageParams(i,:), min(stimDomain));
    averageRandomResponse = makePredicted(model, stimDomain, randomAverageParams(i,:), min(stimDomain));
    
    % Create a new figure for each noise level
    mainFig = figure;
    set(gcf,'Position',[50 50 1200 700]);
    subplot(3,2,[1 3 5]);
    hold on;
    for j = 1:size(qpRows,1)
        params = table2array(qpRows(j,1:length(paramNamesInOrder)));
        qpResponse = makePredicted(model, stimDomain, params, min(stimDomain));
        plot(stimDomain,qpResponse,'-','Color' ,lightQPColor,'LineWidth',.2,'HandleVisibility','off');
    end
    
    for j = 1:size(randRows,1)
        params = table2array(randRows(j,1:length(paramNamesInOrder)));
        randomResponse = makePredicted(model, stimDomain, params, min(stimDomain));
        plot(stimDomain,randomResponse,'-','Color',lightRandColor,'LineWidth',.2,'HandleVisibility','off');
    end
    % Plot veridical and average parameter fits
    plot(stimDomain,averageQPResponse,':','Color',darkQPColor,'LineWidth',4);
    plot(stimDomain,averageRandomResponse,'--','Color',darkRandColor,'LineWidth',4);
    plot(stimDomain,predictedRelativeResponse,'--','Color',veridicalColor,'LineWidth',6);
    %Plot formatting
    ylim([0 1]);
    set(gca,'XScale', 'lin');
    xlabel('Contrast');
    ylabel('Predicted Response, Normalized 0-1');
    legend('Q+','Random','Veridical','Location','Northwest');
    title(horzcat(func2str(model), ' Curves at Noise: ',num2str(noiseLevel), ' SD'));
    hold off;
    
    % Plot histogram for slope
    binRange = linspace(0,1,50);
    hcx = histcounts(qpRows.slope,[binRange Inf]);
    hcy = histcounts(randRows.slope,[binRange Inf]);
    subplot(322);
    b = bar(binRange,[hcx;hcy]','BarWidth',1.5);
    b(1).FaceColor = darkQPColor;
    b(2).FaceColor = darkRandColor;
    hold on;
    plot([sampleSimulatedParams(1) sampleSimulatedParams(1)],ylim,'Color',veridicalColor,'LineWidth',2)
    legend('Q+','Random','Veridical');
    xlabel('Slope Value');
    ylabel('Number of Simulations');
    title('Histogram of Slope Parameter Estimates');
    hold off;
    % Plot histogram for semiSat
    % Plot histogram for slope
    binRange = linspace(0,1,50);
    hcx = histcounts(qpRows.semiSat,[binRange Inf]);
    hcy = histcounts(randRows.semiSat,[binRange Inf]);
    subplot(324);
    hold on;
    b = bar(binRange,[hcx;hcy]','BarWidth',1.5);
    b(1).FaceColor = darkQPColor;
    b(2).FaceColor = darkRandColor;
    plot([sampleSimulatedParams(2) sampleSimulatedParams(2)],ylim,'Color',veridicalColor,'LineWidth',2)
    legend('Q+','Random','Veridical');
    xlabel('Semi-Saturation Value');
    ylabel('Number of Simulations');
    title('Histogram of Semi-Saturation Parameter Estimates');
    hold off;
    


    binRange = linspace(0,2.5,50);
    subplot(326);
    hcx = histcounts(qpRows.maxBOLD,[binRange Inf]);
    hcy = histcounts(randRows.maxBOLD,[binRange Inf]);
    hold on;
    b = bar(binRange,[hcx;hcy]','BarWidth',1.5);
    b(1).FaceColor = darkQPColor;
    b(2).FaceColor = darkRandColor;
    plot([data.maxBOLDSim(1) data.maxBOLDSim(1)],ylim,'Color',veridicalColor,'LineWidth',2)
    legend('Q+','Random','Veridical');
    xlabel('Maximum BOLD Value');
    ylabel('Number of Simulations');
    title(horzcat('Histogram of Maximum BOLD Estimate for Noise = ',num2str(noiseLevel)));
    hold off;
    
    
    set(mainFig,'Units','Inches');
    pos = get(mainFig,'Position');
    set(mainFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(mainFig,horzcat('./',func2str(model),'_',num2str(noiseLevel),'.pdf'),'-dpdf','-r0')

    
end


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

%data.slopeRelToMaxBOLD = (data.slope * data.maxBOLD) ./ data.maxBOLDSim;

SignedResults = varfun(@mean,data,'InputVariables',{'slopeSigned','semiSatSigned','betaSigned','sigmaSigned','maxBOLDSigned'},...
    'GroupingVariables',{'qpPres','sigmaSim'});
 
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

%% Make cool plot
   
function [qpParams,randomParams] = oneNoiseLevel(avgResults,noiseLevel)
% Extract average qpParams and randomParams for one noise level
qpRows = (avgResults.qpPres=="qpControl" & avgResults.sigmaSim==noiseLevel);
a = avgResults(qpRows,:);
qpParams = [a.mean_slope a.mean_semiSat a.mean_beta a.mean_sigma];

randRows = (avgResults.qpPres=="random" & avgResults.sigmaSim==noiseLevel);
b = avgResults(randRows,:);
randomParams = [b.mean_slope b.mean_semiSat b.mean_beta b.mean_sigma];


end

function predictedResponse = makePredicted(model, stimDomain, params, baseline)
predictedResponse = model(stimDomain,params) - model(baseline,params);
end


