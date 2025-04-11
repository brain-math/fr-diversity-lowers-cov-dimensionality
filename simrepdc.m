function [R, E, mu, sigma] = simrepdc(W, targetMean, targetPR, p)
% Sample a pair of R and C with power p
N = size(W, 1);
[R, mu, sigma] = simr(targetMean, targetPR, N);
if length(p) > 1
    k = 1;
else
    k = sqrt(targetMean) * exp(-(1-1/p)*mu-(1-1/p)^2*sigma^2);
end
temp = diag(R.^(1./p-1)) - k * W;
E = svd(temp).^(-2);

end