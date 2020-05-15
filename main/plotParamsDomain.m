function plotParamsDomain(model, paramsDomain, varargin)
%% Plot all possible combos of a parameter domain. 
% TODO: Note - I'm not sure how well this behaves with a 4-parameter model.
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
%   'xParam'                - Integer (Default = 1)
%                             Which parameter to vary along the x subplot axis.
%   'yParam'                - Integer (Default = 2)
%                             Which parameter to vary along the y subplot axis. 
%   'colorParam'            - Integer (Default = 3)
%                             Which parameter to vary color for. 
%   'minStim'               - Scalar (Default = .01)
%                             Lowest value to use for plotting stimulus
%                             domain.
%   'maxStim'               - Scalar (Default = 100)
%                             Highest value to use for plotting stimulus
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

plotParamsDomain(model, paramsDomain);

% Example 2: 
% from early doeSimulate code
paramsDomain = struct;
paramsDomain.Sr = 0.899:0.025:1.099;
paramsDomain.k1 = 0.01:0.04:0.4;
paramsDomain.k2 = 0.01:0.04:0.4;
paramsDomain.beta = 0.8:0.1:1.4; % Amplitude of the scaled response; should converge to unity
paramsDomain.sigma = 0.3:0.2:1;	% Standard deviation of the scaled (0-1) noise


model = @doeTemporalModel;

plotParamsDomain(model, paramsDomain);
%
%}
%% Handle initial inputs
p = inputParser;

% Required input
p.addRequired('model',@(x) isa(x,'function_handle'));
p.addRequired('paramsDomain',@isstruct);


% Optional params
p.addParameter('nOutcomes',51,@isnumeric);
p.addParameter('headroom',.1,@isscalar);
p.addParameter('xParam',1,@isnumeric);
p.addParameter('yParam', 2, @isnumeric);
p.addParameter('colorParam', 3, @isnumeric);
p.addParameter('minStim',.01,@isscalar);
p.addParameter('maxStim',100,@isscalar);
p.addParameter('figWidth',900,@isnumeric);
p.addParameter('figHeight',900,@isnumeric);


% Parse
p.parse( model, paramsDomain, varargin{:});

nOutcomes = p.Results.nOutcomes;
headroom = p.Results.headroom;

% First verify the model is valid and the parameters are accounted for.
[paramNamesInOrder] = checkModel(model,paramsDomain,...,
    'nOutcomes',nOutcomes,'headroom',headroom);

% Models are passed with beta and sigma but we need to ignore them here.
betaIndex = find(strcmp(paramNamesInOrder,'beta'));
sigmaIndex = find(strcmp(paramNamesInOrder,'sigma'));

% Beta and Sigma are not required and won't be plotted for this
nParameters = length(paramNamesInOrder);
% Beta needs to replaced with 1.
paramsDomain.(paramNamesInOrder{betaIndex}) = 1; 

% If sigma is provided, it's just ignored.
if ~isempty(sigmaIndex)
    nParameters = nParameters - 1; 
end

paramSpaceSize = zeros(1,nParameters);
paramVectors = cell(nParameters,1);

for par = 1:nParameters
    paramSpaceSize(par) = length(paramsDomain.(paramNamesInOrder{par}));
    paramVectors{par} = paramsDomain.(paramNamesInOrder{par});
end

% Make one list of all possible parameter combinations
allParameterCombos = allcomb(paramVectors{:});
nAllCombos = length(allParameterCombos);

% Prep for plotting

% Currently 3 or 4 parameters are supported (including Beta)
if nParameters == 4 % including beta
    iterator1 = paramSpaceSize(p.Results.yParam)*paramSpaceSize(p.Results.colorParam);
    iterator2 = paramSpaceSize(p.Results.colorParam);
    colorLength = paramSpaceSize(p.Results.colorParam);
elseif nParameters == 3 % including beta
    iterator1 = paramSpaceSize(p.Results.yParam)*paramSpaceSize(p.Results.colorParam);
    iterator2 = paramSpaceSize(p.Results.colorParam);
    colorLength = paramSpaceSize(p.Results.colorParam);
else
    error('Too many parameters are specified');
end

colorMap = jet(colorLength+1);
color = 1;
subplotNum = 1;
plotRow = 1;
plotColumn = 1;

stimulusFreqHzFine = logspace(log10(p.Results.minStim),log10(p.Results.maxStim),100);

% Initialize figure
fig = figure('Position', [10 10 p.Results.figWidth p.Results.figHeight]);

% Plot all combinations
for c = 1:nAllCombos
    % For the first one, initialize
    if c == 1
        subplot(paramSpaceSize(p.Results.xParam),paramSpaceSize(p.Results.yParam),subplotNum);
        yVals = model(stimulusFreqHzFine,allParameterCombos(c,:));
        semilogx(stimulusFreqHzFine,yVals,'Color',colorMap(color,:));
        plotTitle = sprintf('%s: %.03f',paramNamesInOrder{p.Results.yParam},...,
            paramVectors{p.Results.yParam,1}(plotColumn));
        title(plotTitle); 
        ylabel(sprintf('%s: %.02f',paramNamesInOrder{p.Results.xParam},...,
            paramVectors{p.Results.xParam,1}(plotRow)),'Fontsize',10);
        hold on;

    % every time it iterates through the first parameter (every b*c)
    elseif mod(c-1,iterator1) == 0 
        color = 1;
        plotColumn = 1;
        plotRow = plotRow + 1;
        subplotNum = paramSpaceSize(p.Results.yParam)*(plotRow-1) + plotColumn;
        subplot(paramSpaceSize(p.Results.xParam),paramSpaceSize(p.Results.yParam),subplotNum);
        ylabel(sprintf('%s: %.02f',paramNamesInOrder{p.Results.xParam},...
            paramVectors{p.Results.xParam,1}(plotRow)),'Fontsize',10);
        hold on;
    elseif mod(c-1,iterator2) == 0
        color = 1;
        plotColumn = plotColumn + 1;
        subplotNum = paramSpaceSize(p.Results.yParam)*(plotRow-1) + plotColumn;
        subplot(paramSpaceSize(p.Results.xParam,1),paramSpaceSize(p.Results.yParam),subplotNum);
        hold on;
    end
    
    if plotRow == 1
        plotTitle = sprintf('%s: %.03f',paramNamesInOrder{p.Results.yParam},...,
            paramVectors{p.Results.yParam,1}(plotColumn));
        title(plotTitle); 
    end
    
    color = color + 1;
    yVals = model(stimulusFreqHzFine,allParameterCombos(c,:));
    semilogx(stimulusFreqHzFine,yVals,'Color',colorMap(color,:));
    set(gca,'xtick',[],'ytick',[],'XScale', 'log')
    
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';


colormap(jet(colorLength+1));
cbh = colorbar('YTickLabel',[0 round(paramVectors{p.Results.colorParam},2)]);
colorTitleHandle = get(cbh,'Title');
titleString = sprintf('%s',paramNamesInOrder{p.Results.colorParam});
set(colorTitleHandle ,'String',titleString,'Fontsize',15);
set(cbh, 'Position', [.07 .2 .02 .6]);
hold off;
    