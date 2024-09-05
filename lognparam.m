function [mu, sigma] = lognparam(targetMean, targetPR)
%LOGNPARAM Determine the appropriate parameters of the log normal 
% distribution for the desired mean and PR values
sigmaSq = -log(targetPR);
sigma = sqrt(sigmaSq);
mu = log(targetMean) - sigmaSq/2;

end
