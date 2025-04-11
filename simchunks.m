function [r, r0, prRate, prCov, prCov0, prV, R, E, C] = simchunks(nChunk, W, targetMean, prRange)
% Simulate chunks of data (i.e. different input patterns)
N = size(W, 1); N0 = N / 2;
prRate = zeros(nChunk, 1); prCov = prRate; prV = prRate;
prRate0 = zeros(nChunk, 2); prCov0 = prRate0;

for iChunk = 1 : nChunk
    targetPR = rand*(prRange(2)-prRange(1)) + prRange(1);
    R = simr(targetMean, targetPR, N);
    temp = diag(R.^(-1/2)) - W; E = svd(temp) .^ (-2);
    temp = inv(temp); C = temp * temp';
    prRate(iChunk) = pr(R);
    prCov(iChunk) = pr(E);
    prV(iChunk) = pr(diag(C));
    prRate0(iChunk, 1) = pr(R(1:N0)); 
    prCov0(iChunk, 1) = pr(eig(C(1:N0, 1:N0)));
    prRate0(iChunk, 2) = pr(R(1+N0:end)); 
    prCov0(iChunk, 2) = pr(eig(C(1+N0:end, 1+N0:end)));
end

r = corrcoef(prRate, prCov); r = r(1, 2);
r0 = zeros(1, 2);
temp = corrcoef(prRate0(:, 1), prCov0(:, 1)); r0(1) = temp(1, 2);
temp = corrcoef(prRate0(:, 2), prCov0(:, 2)); r0(2) = temp(1, 2);

end

