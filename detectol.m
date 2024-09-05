function [nOl, eigSim] = detectol(E, R, T, fitted, nRepeat, p)
%DETECTOL Detect outliers for one covariance matrix
S = length(R);
W = randn(fitted.N) * fitted.g / sqrt(fitted.N);
Z = randn(fitted.N, T) / sqrt(T);
eigSim = fitted.k * repeatedsample(W, Z, R, nRepeat);
p0 = 1 - (1 - p)^(1/S);
nOl  = markol(E, eigSim, p0);

end

function nOl = markol(E, eigSim, p)
% Mark the outliers from the fitted model
eigSim = sort(eigSim);
l = length(eigSim);
pQuantile = eigSim(floor((1-p)*l));
nOl = nnz(E > pQuantile);

end
