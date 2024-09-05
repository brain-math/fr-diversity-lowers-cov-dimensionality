% Data analysis pipeline for Tian et al. (2024)
% Code by: Gengshuo John Tian
% Requires public V1-V2 data from Zandvakili & Kohn (2015)
% V1-V2 data download: https://crcns.org/data-sets/vc/v1v2-1/about_v1v2-1

set(groot, 'defaultLineLineWidth', 2)
set(groot, 'defaultAxesLineWidth', 2)
set(groot, 'defaultAxesFontSize', 20)
set(groot, 'defaultAxesTickDir', 'out')
set(groot,  'defaultAxesTickDirMode', 'manual')

% Put downloaded data in the same directory
addpath ./v1v2/v1-v2_gratings/software/
dataDir = './v1v2/v1-v2_gratings/mat_neural_data/'; 

%% Preprocess example session
param.nCond = 8; % 8 stimulus conditions (drifting grating orientations)
param.binSize = 100; % Spike count bin size.
param.winLength = 1000; % The length of time window used in each trial

load([dataDir, '106r002p70.mat'])
areaData = procsession(neuralData, param);

%% Fitting a null model to the eigenvalue distribution to identify outliers
% fitParam.p = 0.05; % p-value for detecting outliers
% fitParam.nRepeat = 1e2; % See fitkng.m for definition
% fitParam.verbFlag = 1; % Verbose output during fitting
% fitParam.plotFlag = 1; % Plot the Bayesian optimization process
% fitParam.maxIter = 150; % Maximum iterations for fitting
% fitParam.maxN = 2000; % Maximum underlying population size
% fitted = fitarea(areaData, fitParam); % Fit model parameters k, N, g
% nOl = getol(areaData, fitted, fitParam); % Detect outliers

% The above code for fitting a null model is commented out because it takes
% considerable time to run. We provide here the result of the fitting (the
% number of eigenvalue outliers) directly. To run the optimization
% yourself to get nOl, simply uncomment the above code.
nOl = 3;

%% Calculate and plot D(C) and D(R) (Fig. 3b, left)
[DR, DC] = getpr(areaData, nOl);

figure; plot(DR, DC, 'o')
xlabel('D(R)'); ylabel('D(C)')
