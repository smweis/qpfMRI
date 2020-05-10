function [] = plotParamsDomain(model, params)
% model     -          functionHandle.
%
% params    -          A struct consisting of all parameter domains. 
%                      Current functionality can plot 3 parameters.
%                      DoE    (n=5): Sr, k1, k2, beta, sigma
%                      Watson (n=5): tau, kappa, zeta, beta, sigma
%
% Example: 
%{
% Params domain to plot
params = struct;
params.Sr = logspace(log10(0.01),log10(5),10);
params.k1 = linspace(.01,.08,10);
params.k2 = logspace(log10(0.01),log10(100),10);

params.beta = 1;

model = @doeTemporalModel;

plotParamsDomain(model, params);
%
%}

% Model specific processing. This should be functionalized. 
modelsCreated = {'doeTemporalModel','watsonTemporalModel'};
s = functions(model);
assert(any(strcmp(modelsCreated,s.function)==1),'Model not defined');

if contains(s.function,'doe')
    assert(isfield(params,'Sr'),'Sr missing or misnamed for DoE model.');
    assert(isfield(params,'k1'),'k1 missing or misnamed for DoE model.');
    assert(isfield(params,'k2'),'k2 missing or misnamed for DoE model.');
elseif contains(s.function,'watson')
    assert(isfield(params,'tau'),'tau missing or misnamed for Watson model.');
    assert(isfield(params,'kappa'),'kappa missing or misnamed for Watson model.');
    assert(isfield(params,'zeta'),'zeta missing or misnamed for Watson model.');
end

parameterNames = fieldnames(params);
nParameters = length(parameterNames);
paramSpaceSize = zeros(1,nParameters);
paramVectors = cell(nParameters,1);

for p = 1:nParameters
    paramSpaceSize(p) = length(params.(parameterNames{p}));
    paramVectors{p} = params.(parameterNames{p});
end

% Make one list of all possible parameter combinations
allParameterCombos = allcomb(paramVectors{:});
nAllCombos = length(allParameterCombos);

% Prep for plotting


if nParameters == 4 % including beta
    iterator1 = paramSpaceSize(2)*paramSpaceSize(3);
    iterator2 = paramSpaceSize(3);
    colorLength = paramSpaceSize(3);
elseif nParameters == 3 % including beta
    iterator1 = paramSpaceSize(1)*paramSpaceSize(2);
    iterator2 = paramSpaceSize(2);
    colorLength = paramSpaceSize(2);
end

colorMap = jet(colorLength+1);
color = 1;
subplotNum = 1;
plotRow = 1;
plotColumn = 1;
stimulusFreqHzFine = logspace(log10(0.01),log10(100),100);

% Initialize figure
fig = figure;

% Plot all combinations
for c = 1:nAllCombos
    % For the first one, initialize
    if c == 1
        subplot(paramSpaceSize(1),paramSpaceSize(2),subplotNum);
        yVals = model(stimulusFreqHzFine,allParameterCombos(c,:));
        semilogx(stimulusFreqHzFine,yVals,'Color',colorMap(color,:));
        plotTitle = sprintf('%s: %.03f',parameterNames{2},paramVectors{2,1}(plotColumn));
        title(plotTitle); 
        ylabel(sprintf('%s: %.02f',parameterNames{1},paramVectors{1,1}(plotRow)),'Fontsize',10);
        hold on;

    % every time it iterates through the first parameter (every b*c)
    elseif mod(c-1,iterator1) == 0 
        color = 1;
        plotColumn = 1;
        plotRow = plotRow + 1;
        subplotNum = paramSpaceSize(2)*(plotRow-1) + plotColumn;
        subplot(paramSpaceSize(1),paramSpaceSize(2),subplotNum);
        ylabel(sprintf('%s: %.02f',parameterNames{1},paramVectors{1,1}(plotRow)),'Fontsize',10);
        hold on;
    elseif mod(c-1,iterator2) == 0
        color = 1;
        plotColumn = plotColumn + 1;
        subplotNum = paramSpaceSize(2)*(plotRow-1) + plotColumn;
        subplot(paramSpaceSize(1),paramSpaceSize(2),subplotNum);
        hold on;
    end
    
    if plotRow == 1
        plotTitle = sprintf('%s: %.03f',parameterNames{2},paramVectors{2,1}(plotColumn));
        title(plotTitle); 
    end
    
    color = color + 1;
    yVals = model(stimulusFreqHzFine,allParameterCombos(c,:));
    semilogx(stimulusFreqHzFine,yVals,'Color',colorMap(color,:));
    set(gca,'xtick',[],'ytick',[])
    
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';


colormap(jet(colorLength+1));
cbh = colorbar('YTickLabel',[0 round(paramVectors{end-1},2)]);
colorTitleHandle = get(cbh,'Title');
titleString = sprintf('%s',parameterNames{end-1});
set(colorTitleHandle ,'String',titleString,'Fontsize',15);
set(cbh, 'Position', [.07 .2 .02 .6]);
hold off;
    
    