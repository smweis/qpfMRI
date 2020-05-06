function [outputFilename] = genSimCsv(params, paramToVary, numSims)
%% Creates csv containing all paramater combinations to use for a batch of simulations
%   Input:
%       params      - struct containing all assumed parameter values
%       paramToVary - string specifying which parameter will be varied in
%                     the csv
%       numSims     - integer value representing the number of simulations each
%                     value of paramToVary will be ran
%
%   Output:
%       outputFileName: string containing the absolute path to the
%                       generated csv file
%  
%   Example: genSimCsv(samplePrams, "k1", 5)
%
%   General thoughts:
%     1. Only the variable we are evaluating should have a domain specified in $params 
%     2. Should the entire domain of this variable be represented in the csv?
%     
numParams = numel(fieldnames(params))-1;

for i = params.(paramToVary+"Domain")
    for j = 1:numSims
        dlmwrite("csvtest.csv", i, '-append');
    end
end
end

