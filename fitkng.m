function [k, N, g, cost, results] = fitkng(E, R, T, fitParam)
%FITKNG Fit a model to the eigenvalues
%   Arguments:
%       E: covariance matrix eigenvalues from data to be fitted to.
%       R: vector of firing rates.
%       T: number of temporal samples.
%       fitParam.nRepeat: sample from the big covariance matrix this many 
%               times to get a good estimate of the sampled covariance 
%               matrix's eigenvalue distribution.
%   Returns:
%       k: scaling factor.
%       N: underlying population size.
%       g: variance parameter of the connection matrix.
%       cost: final cost function value.

meanR = cellfun(@mean, R); lengthR = cellfun(@length, R);
meanE = mean(cellfun(@mean, E));
lb = [ceil(1.1*max(lengthR)), 0]; ub = [fitParam.maxN, 0.9/sqrt(max(meanR))];
meanR = mean(meanR);
if fitParam.plotFlag
    plotFcn = {@plotObjectiveModel, @plotMinObjective};
else
    plotFcn = [];
end

% Bayesian optimization
varN = optimizableVariable('N', [lb(1), ub(1)], 'Type', 'integer');
varG = optimizableVariable('g', [lb(2), ub(2)]);
fun = @(x) obj(E, R, T, x.N, x.g, meanR, meanE, fitParam.nRepeat);
results = bayesopt(fun, [varN, varG], 'Verbose', fitParam.verbFlag,...
    'PlotFcn', plotFcn, 'MaxObjectiveEvaluations', fitParam.maxIter);

N = results.XAtMinEstimatedObjective.N; 
g = results.XAtMinEstimatedObjective.g;
mu = meanR / (1 - g^2*meanR); k = meanE / mu;
cost = results.MinEstimatedObjective;

end

function f = obj(E, R, T, N, g, meanR, meanE, nRepeat)
W = randn(N) * g / sqrt(N);
nChunk = length(E);
mu = meanR / (1 - g^2*meanR); k = meanE / mu;

f = 0;
for iChunk = 1 : nChunk
    Z = randn(N, T{iChunk}) / sqrt(T{iChunk});
    sim = repeatedsample(W, Z, R{iChunk}, nRepeat);
    sim = k * sim;
    f = f + ad(E{iChunk}, sim);
end

end

