function [fig] = plotParamsDomain(myQpfmriParams, varargin)
%% Plot all possible combos of a parameter domain. 
% Plots all possible combinations of a parameter domain across a range of
% stimulus values.
% Inputs:
%   myQpfmriParams          - Struct. Set of parameters used for qpfmri.
%                             See qpfmriParams function for more details
% Optional key/value pairs:
%   'xParam'                - Integer (Default = 1)
%                             Which parameter to vary along the x subplot axis.
%   'yParam'                - Integer (Default = 2)
%                             Which parameter to vary along the y subplot axis. 
%   'colorParam'            - Integer (Default = 3)
%                             Which parameter to vary color for. 
%                             domain.
%   'figWidth'              - Integer (Default = 900)
%                             Width of figure window size.
%   'figHeight'             - Integer (Default = 900)
%                             Height of figure window size.
%
% Example: 
%{
% Params domain to plot
% Example 1:
params = struct;
params.Sr = logspace(log10(0.01),log10(5),10);
params.k1 = linspace(.01,.08,10);
params.k2 = logspace(log10(0.01),log10(100),10);
paramsDomain.beta = 0.8:0.1:1.4; % Amplitude of the scaled response; should converge to unity
paramsDomain.sigma = 0.3:0.2:1;	% Standard deviation of the scaled (0-1) noise
model = @doeTemporalModel;
stimulusDomain = [1.875, 3.75, 7.5, 15, 30, 60];

plotParamsDomain(myQpfmriParams);

% Example 2: 
% from early doeSimulate code
paramsDomain = struct;
paramsDomain.Sr = 0.899:0.025:1.099;
paramsDomain.k1 = 0.01:0.04:0.4;
paramsDomain.k2 = 0.01:0.04:0.4;
paramsDomain.beta = 0.8:0.1:1.4; % Amplitude of the scaled response; should converge to unity
paramsDomain.sigma = 0.3:0.2:1;	% Standard deviation of the scaled (0-1) noise
model = @doeTemporalModel;
stimulusDomain = logspace(log10(0.01),log10(200),100);;
stimDomainSpacing = 'log';

myQpfmriParams = (model,paramsDomain,'stimulusDomain',stimulusDomain,...,
'stimulusDomainSpacing',stimulusDomainSpacing);

plotParamsDomain(myQpfmriParams);

%Example 3

model = @logistic;

paramsDomain = struct;
paramsDomain.slope = linspace(.01,1,40);
paramsDomain.semiSat = linspace(.01,1,40);
paramsDomain.beta = 0.8:0.1:1.4; 
paramsDomain.sigma = linspace(.3,1.5,8);

stimulusDomain = {linspace(.01,1,10)};
stimulusDomainSpacing = 'lin';

myQpfmriParams = (model,paramsDomain,'stimulusDomain',stimulusDomain,...,
'stimulusDomainSpacing',stimulusDomainSpacing);

plotParamsDomain(myQpfmriParams);

%}
%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('myQpfmriParams',@isstruct);

% Optional params
p.addParameter('xParam',1,@isnumeric);
p.addParameter('yParam', 2, @isnumeric);
p.addParameter('colorParam', 3, @isnumeric);
p.addParameter('figWidth',900,@isnumeric);
p.addParameter('figHeight',900,@isnumeric);


% Parse
p.parse( myQpfmriParams, varargin{:});

% First verify the model is valid and the parameters are accounted for.
[myQpfmriParams] = checkModel(myQpfmriParams);


% Beta and Sigma are not required and won't be plotted for this
nParameters = length(myQpfmriParams.paramNamesInOrder);
% Beta needs to replaced with 1.
myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{betaIndex}) = 1; 

% If sigma is provided, it's just ignored.
if ~isempty(myQpfmriParams.sigmaIndex)
    nParameters = nParameters - 1; 
end

paramSpaceSize = zeros(1,nParameters);
paramVectors = cell(nParameters,1);

for par = 1:nParameters
    
    paramVectors{par} = myQpfmriParams.paramsDomain.(myQpfmriParams.paramNamesInOrder{par});
    if length(paramVectors{par}) > 10
        warning('Too many divisions for %s to plot in a reasonable way. Parameter domain made coarser.',...,
            myQpfmriParams.paramNamesInOrder{par});
        paramVectors{par} = linspace(min(paramVectors{par}),max(paramVectors{par}),10);
    end
    paramSpaceSize(par) = length(paramVectors{par});
end


% Currently 3 or 4 parameters are supported (including Beta)
if nParameters == 4  % including beta
    figWidth = p.Results.figWidth;
    figHeight = p.Results.figHeight;
    colorLength = paramSpaceSize(p.Results.colorParam);
    iterator1 = paramSpaceSize(p.Results.yParam)*paramSpaceSize(p.Results.colorParam);
    iterator2 = paramSpaceSize(p.Results.colorParam);
    xParam = paramSpaceSize(p.Results.xParam);
    yParam = paramSpaceSize(p.Results.yParam);
    plotTitleParam = p.Results.yParam;
elseif nParameters == 3
    if p.Results.colorParam == 3
        colorParam = p.Results.yParam;
    end
    figWidth = 500;
    figHeight = 1000;
    xParam = paramSpaceSize(p.Results.xParam);
    yParam = 1;
    colorLength = paramSpaceSize(p.Results.yParam);
    iterator1 = paramSpaceSize(p.Results.xParam)*paramSpaceSize(p.Results.yParam);
    iterator2 = paramSpaceSize(p.Results.yParam);
    plotTitleParam = p.Results.xParam;
    
else
    error('Too many parameters are specified');
end




% Make one list of all possible parameter combinations
allParameterCombos = allcomb(paramVectors{:});
nAllCombos = length(allParameterCombos);

% Prep for plotting

if strcmpi(myQpfmriParams.stimulusDomainSpacing,'log')
    plotFunc = @semilogx;
else 
    plotFunc = @plot;
end


% Initialize figure
fig = figure('Position', [10 10 figWidth figHeight]);


colorMap = jet(colorLength+1);
color = 1;
subplotNum = 1;
plotRow = 1;
plotColumn = 1;

% Plot all combinations
for c = 1:nAllCombos
    % For the first one, initialize
    if c == 1
        subplot(xParam,yParam,subplotNum);
        yVals = myQpfmriParams.model(myQpfmriParams.stimulusDomain{:},allParameterCombos(c,:));
        plotFunc(myQpfmriParams.stimulusDomain{:},yVals,'Color',colorMap(color,:));
        plotTitle = sprintf('%s: %.03f',paramNamesInOrder{plotTitleParam},...,
            paramVectors{plotTitleParam,1}(plotColumn));
        title(plotTitle); 
        hold on;

    % every time it iterates through the first parameter (every b*c)
    elseif mod(c-1,iterator1) == 0 
        color = 1;
        plotColumn = 1;
        plotRow = plotRow + 1;
        subplotNum = paramSpaceSize(p.Results.yParam)*(plotRow-1) + plotColumn;
        subplot(xParam,yParam,subplotNum);
        hold on;
    elseif mod(c-1,iterator2) == 0
        color = 1;
        plotColumn = plotColumn + 1;
        
        if nParameters == 4
            ylabel(sprintf('%s: %.02f',paramNamesInOrder{p.Results.xParam},...
                paramVectors{p.Results.xParam,1}(plotRow)),'Fontsize',10);
        end
            
        subplotNum = paramSpaceSize(p.Results.yParam)*(plotRow-1) + plotColumn;
        subplot(xParam,yParam,subplotNum);
        hold on;
    end
    
    if plotRow == 1
        plotTitle = sprintf('%s: %.03f',paramNamesInOrder{plotTitleParam},...,
            paramVectors{plotTitleParam,1}(plotColumn));
        title(plotTitle); 
    end

    color = color + 1;
    yVals = myQpfmriParams.model(myQpfmriParams.stimulusDomain{:},allParameterCombos(c,:));
    plotFunc(myQpfmriParams.stimulusDomain{:},yVals,'Color',colorMap(color,:));
    
    if strcmpi(myQpfmriParams.stimulusDomainSpacing,'log')
        set(gca,'xtick',[],'ytick',[],'XScale', 'log');
    end
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';


colormap(jet(colorLength+1));
cbh = colorbar('YTickLabel',[0 round(paramVectors{colorParam},2)]);
colorTitleHandle = get(cbh,'Title');
titleString = sprintf('%s',paramNamesInOrder{colorParam});
set(colorTitleHandle ,'String',titleString,'Fontsize',15);
set(cbh, 'Position', [.07 .2 .02 .6]);
drawnow;
hold off;
    