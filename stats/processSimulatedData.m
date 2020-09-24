function [stats] = processSimulatedData(resultsDir,groupingVars)


posterFormat = struct;
posterFormat.lightRandColor = '#3D92C9';
posterFormat.darkRandColor = '#165172';
posterFormat.lightQPColor = '#F99C16';
posterFormat.darkQPColor = '#FA4515';
posterFormat.veridicalColor = '#000000';

% Grab the names of all the mat files in the results directory
resultsNames = dir(fullfile(resultsDir,'*.mat'));

% Check if data have already been processed.
if ~isempty(find(strcmp({resultsNames.name}, 'loadedData.mat')==1,1))
    load(fullfile(resultsDir,'loadedData.mat'));
% If the number of mat files in the results directory doesn't match the
% number of files found in loadedData.mat, re-run it to grab the right
% data.
    if size(data.qpfmriResults,1) ~= size(resultsNames,1)-1
        data = loadSimulatedData(resultsDir);
    end
else 
    data = loadSimulatedData(resultsDir);
end

% Turn the struct into a table for easier plotting.
dataTable = struct2table(data.qpfmriResults);

% Grab one of files. Useful for a few details that aren't likely to change
% across simulations. 
myQpfmriParams = load(fullfile(resultsNames(2).folder,resultsNames(2).name));
myQpfmriParams = myQpfmriParams.qpfmriResults;

% Get rid of data, which is enormous. 
clear('data');

% Re-code a couple variables that don't store properly.
for i = 1:height(dataTable)
    temp = dataTable.entropyOverTrials(i,end);
    dataTable.entropyFinal(i) = temp{end}(end);
    dataTable.maxBOLDFinalPreFit(i) = dataTable.maxBOLDoverTrials(i,end);
    if contains(dataTable.outNum{i},'qp')
        dataTable.qpPresFinal{i} = 'qp';
    else
        dataTable.qpPresFinal{i} = 'random';
    end
end



inputVars = {'nOutcomes','noiseSD','qpPresFinal','maxBOLDFinal',...,
    'maxBOLDFinalPreFit','psiParamsBadsFinal','entropyFinal'};

dsa = dataTable(:,inputVars);
G = findgroups(dataTable.qpPresFinal,dataTable.nOutcomes,dataTable.noiseSD);
stats = grpstats(dsa,groupingVars,{'mean','median','std','meanci'});
qpRows = (dataTable.qpPresFinal=="qp" & dataTable.noiseSD==.1);
qpRows = dataTable(qpRows==1,:);
randRows = (dataTable.qpPresFinal=="random" & dataTable.noiseSD==.1);
randRows = dataTable(randRows==1,:);

modelResponseNoCorrection = myQpfmriParams.model(myQpfmriParams.stimulusDomain{:},myQpfmriParams.simulatedPsiParams);




% Create a new figure for each level
mainFig = figure('Position', get(0, 'Screensize'));
set(gcf,'color','w');


hold on;

for j = 1:size(randRows,1)
    params = randRows.psiParamsBadsFinal(j,:);
    randomResponse = makePredicted(myQpfmriParams.model, myQpfmriParams.stimulusDomain{:}, params, min(myQpfmriParams.stimulusDomain{:}),randRows.maxBOLDFinal(j));
    simResponse = myQpfmriParams.model(myQpfmriParams.stimulusDomain{:},params);
    randRows.rSquared(j) = corr(modelResponseNoCorrection',simResponse')^2;
    plot2 = plot(myQpfmriParams.stimulusDomain{:},randomResponse,'-','Color',posterFormat.darkRandColor,'LineWidth',2,'HandleVisibility','off');
    plot2.Color(4) = 0.1;
end

for j = 1:size(qpRows,1)
    params = qpRows.psiParamsBadsFinal(j,:);
    qpResponse = makePredicted(myQpfmriParams.model, myQpfmriParams.stimulusDomain{:}, params, min(myQpfmriParams.stimulusDomain{:}),qpRows.maxBOLDFinal(j));
    simResponse = myQpfmriParams.model(myQpfmriParams.stimulusDomain{:},params);
    qpRows.rSquared(j) = corr(modelResponseNoCorrection',simResponse')^2;
    plot2 = plot(myQpfmriParams.stimulusDomain{:},qpResponse,'-','Color',posterFormat.darkQPColor,'LineWidth',2,'HandleVisibility','off');
    plot2.Color(4) = 0.1;
end

a = 'hi';

function predictedResponse = makePredicted(model, stimDomain, params, baseline, maxBOLD)
predictedResponse = model(stimDomain,params) - model(baseline,params);
%predictedResponse = predictedResponse .* maxBOLD;
%predictedResponse = model(stimDomain,params);
end

end
