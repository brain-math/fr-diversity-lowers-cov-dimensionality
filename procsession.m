function areaData = procsession(neuralData, param)
%PROCSESSION Process a session of data
[spikesDriven, stim] = ExtractSpikes(neuralData, param.binSize, 'TrialPeriod', 'Driven');
spikesDriven = spikesDriven{1}; % Use V1 data
spikesBlank = ExtractSpikes(neuralData, param.binSize, 'TrialPeriod', 'Blank');
spikesBlank = spikesBlank{1}; % Use V1 data
nTrial = size(spikesBlank, 3); trialPerChunk = floor(nTrial/param.nCond);
nBin = param.winLength / param.binSize; % Number of bins per trial

n = 1;
for iCond = 1 : param.nCond
    X = spikesDriven(:, end-nBin+1:end, stim==iCond);
    areaData(n) = procx(X, param); n = n + 1;
    % Parse blank screen data according to time instead of preceding stimulus
    X = spikesBlank(:, end-nBin+1:end, (iCond-1)*trialPerChunk+1:iCond*trialPerChunk);
    areaData(n) = procx(X, param); n = n + 1;
end

end

function temp = procx(X, param)
% Process spike count data X.
X = reshape(X, [size(X, 1), size(X, 2)*size(X, 3)])';
meanX = mean(X); Xcentered = X - meanX; R = meanX / param.binSize * 1e3;

% Calculate eigenvalues of C using svd
S = svd(Xcentered);
E = S.^2 / size(Xcentered, 1) / param.binSize * 1e3;
temp.X = X; temp.R = R'; temp.E = E;

end