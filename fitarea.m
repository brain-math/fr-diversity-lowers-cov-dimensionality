function [fitted, cost, results] = fitarea(areaData, fitParam)
%FITAREA Fit g and N for one areaData structure
[E, R, T] = gatherer(areaData);
[k, N, g, cost, results] = fitkng(E, R, T, fitParam);
fitted.k = k; fitted.N = N; fitted.g = g;

end