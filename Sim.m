% Numerical simulations for Tian et al. (2024)
% Code by: Gengshuo John Tian

set(groot, 'defaultLineLineWidth', 2)
set(groot, 'defaultAxesLineWidth', 2)
set(groot, 'defaultAxesFontSize', 20)
set(groot, 'defaultAxesTickDir', 'out')
set(groot,  'defaultAxesTickDirMode', 'manual')
co = colororder;

rng(1) % Set random seed

%% Main result (Fig. 2a)
N = 1e3; % Number of neurons (5e3 in the paper)
nR = 200; % Number of conditions (i.e. operating points)
targetMean = 10; % Population mean firing rate
g = 0.15; % Weight distribution variance parameter
DC = zeros(1, nR); DR = zeros(1, nR); % Initialize D(C) and D(R)

h = waitbar(0);
W = randn(N) * g / sqrt(N); % Generate weight matrix W
for iR = 1 : nR
    targetPR = rand * 0.8 + 0.1; % Sample D(R) uniformly between 0.1 and 0.9
    [R, E] = simrc(W, targetMean, targetPR);
    DR(iR) = pr(R); DC(iR) = pr(E);
    waitbar(iR/nR, h)
end
delete(h)

x = 0 : 0.1 : 1; y = x * (1 - g^2 * targetMean)^2; % Theoretical prediction

% Plot result
figure; plot(DR, DC, 'o'); hold on
plot(x, y); legend('Simulated data', 'Theory')
xlabel('D(R)'); ylabel('D(C)')

%% Two-population model (Figs. 1d, S1a)
gRatio = 3; % g22 / g11
N0 = 500; % Number of neurons in each population
N = N0 * 2; % Total number of neurons
targetMean = 10; % Population mean firing rate
gm = 1 / sqrt(targetMean); % Maximum g allowed in the linearizable regime

g0 = 0.6 * gm; % Total connectivity variance parameter
DR = 0.6; % D(R)
% Scan through the above two parameters to get Fig. S1a

% Calculate the g for each block
g2total = g0^2 * 4;
g1 = sqrt(g2total / (gRatio^2 + 3)); % g11 = g12 = g21 = g1
g2 = g1 * gRatio; % g22 = g2

DC0 = DR * (1 - g0^2 * targetMean) ^ 2; % Theoretical D(C)
W = [randn(N0)*g1, randn(N0)*g1; randn(N0)*g1, randn(N0)*g2] / sqrt(N);
R = simr(targetMean, DR, N); % Generate firing rates
temp = diag(R.^(-1/2)) - W; E = svd(temp) .^ (-2); % Get covariance spectrum
DC = pr(E); % D(C) from simulation. Compare with DC0

% Firing rate distribution
figure; histogram(R); xlabel('Firing rate (Hz)')

% Visualize the correlation matrix
C = inv(temp); C = C * C'; % Calculate the covariance matrix
D = sqrt(diag(C)); Cor = C ./ (D * D');  % Calculate the correlation matrix
Cor(logical(eye(N))) = mean(triu_ele(Cor)); % Suppress the diagonal
figure; imagesc(Cor); title('Correlation matrix')
cm = max(abs(Cor(:))); clim([-cm, cm]); colormap gray; colorbar

% Eigenvalue distribution
x = linspace(0.1, 300, 1e4);
p = computepdf(x, R, g0); 
c = cumsum(p); c = c / c(end);
nRank = 100; Edata = sort(E, 'descend'); Edata = Edata(1 : nRank);
ratio = 10; nRankT = nRank * ratio; Etheory = zeros(1, nRankT);
for iRank = 1 : nRankT
    Etheory(iRank) = x(find(c>1-iRank/ratio/N, 1));
end
figure; plot(1/ratio:1/ratio:nRank, Etheory); hold on
scatter(1:nRank, Edata, 'filled'); xlabel('Rank'); ylabel('Eigenvalue')

%% [D(R), D(C)] correlation in the two-population model (Figs. S1b)
nRepeat = 100; % Repeat to generate histogram
nChunk = 20; % Number of conditions
rAll = zeros(nRepeat, 1); % [D(R), D(C)] correlation for both populations
r0All = zeros(nRepeat, 1); % [D(R), D(C)] correlation for each population
N0 = 100; % Number of neurons in each population
N = N0 * 2; % Total number of neurons
prRange = [0.3, 0.7]; % A similar range of D(R) as in real data
targetMean = 10; % Population mean firing rate
gm = 1 / sqrt(targetMean); % Maximum g allowed in the linearizable regime
g0 = 0.6 * gm; % Total connectivity variance parameter

ratio1 = 3; % g22 / g11
ratio2 = 2; % g11 / g12
% Scan through the above two parameters to get Fig. S1b

% Calculate the g for each block
g2total = g0^2 * 4;
g12 = sqrt(g2total / (2+ratio2^2+(ratio1*ratio2)^2)); % g21 = g12
g11 = g12 * ratio2; 
g22 = g11 * ratio1;

h = waitbar(0); 
for iRepeat = 1 : nRepeat
    W = [randn(N0)*g11, randn(N0)*g12; randn(N0)*g12, randn(N0)*g22] / sqrt(N);
    [r, r0] = simchunks(nChunk, W, targetMean, prRange);
    rAll(iRepeat) = r; 
    r0All(iRepeat) = mean(r0);
    waitbar(iRepeat/nRepeat, h)
end
delete(h); 

x = -1 : 0.1 : 1;
figure; histogram(r0All, x); hold on; histogram(rAll, x)
xlabel('Corr. of D(R) and D(C)'); 
legend('One population', 'Two populations')

%% Spatial sampling (Fig. S2)
N = 1e3; g = 0.2; targetMean = 10; nR = 100;
Sall = [48, 200, 500, 1000]; % Subpopulation size
nS = length(Sall);
W = randn(N) * g / sqrt(N);

x = 0.01 : 0.01 : 0.99; nx = length(x); % For plotting theory
h = waitbar(0); t = 0; lStr = {}; lH = []; figure
for iS = 1 : nS
    S = Sall(iS);
    sampled = false(1, N);
    sampled(1:S) = true;
    prRate = zeros(nR, 1); 
    prCov = zeros(nR, 1);
    for iR = 1 : nR
        targetPR = rand * 0.8 + 0.1;
        R = simr(targetMean, targetPR, N);
        temp = inv(diag(R.^(-1/2)) - W);
        C = temp * temp';
        E = eig(C(sampled, sampled));
        prRate(iR) = pr(R(sampled));
        prCov(iR) = pr(E);
        t = t + 1; waitbar(t/nR/nS, h)
    end
    lH(iS) = scatter(prRate, prCov, [], co(iS, :), 'filled'); hold on
    if iS == nS % No sampling case
        prCovPred = x * (1 - g^2 * targetMean)^2;
        lStr(iS) = 'f=\infty';
    else
        f = S / (N-S); g2sq = g^2/(1+f); g1sq = f*g2sq;
        factor = 1/(1/(1-f*g2sq*targetMean)^2+f*(1/(1-g2sq*targetMean)-1)^2);
        prCovPred = factor * x;
        lStr(iS) = ['f=', num2str(f)];
    end
    plot(x, prCovPred, 'color', co(iS, :))
end
delete(h)
xlabel('D(R)'); ylabel('D(C)'); legend(lH, lStr)

%% Different power p of the transfer function (Fig. S3a)
N = 1e3; targetMean = 10; g = 0.5 / sqrt(targetMean);
P = [0.5, 1, 2, 6]; % Scan through this to get Fig. S3b
nP = length(P); prRange = [0.1, 0.9]; nR = 100;
x = 0.01 : 0.01 : 0.99; nx = length(x);

t = 0; h = waitbar(0); lH = []; lStr = {}; figure
for iP = 1 : nP
    p = P(iP);
    W = randn(N) * g / sqrt(N);
    prRate = zeros(nR, 1); prCov = prRate;
    for iR = 1 : nR
        targetPR = rand*(prRange(2)-prRange(1))+prRange(1);
        [R, E] = simrepdc(W, targetMean, targetPR, p);
        prRate(iR) = pr(R); prCov(iR) = pr(E);
        t = t + 1; waitbar(t/nP/nR, h)
    end
    lH(iP) = scatter(prRate, prCov, [], co(iP, :), 'filled'); hold on
    lStr{iP} = ['p=', num2str(p)];

    prPred = zeros(nx, 1);
    for ix = 1 : nx
        targetPR = x(ix);
        [~, sigma] = lognparam(targetMean, targetPR);
        prPred(ix) = exp(-4*(p-1)^2*sigma^2/p^2) * (1 - g^2 * targetMean)^2;
    end
    plot(x, prPred, 'color', co(iP, :))
end
delete(h)
xlabel('D(R)'); ylabel('D(C)'); legend(lH, lStr)

%% A distribution of power p across neurons (Fig. S3c)
N = 1e3; targetMean = 10; g = 0.5 / sqrt(targetMean);
S = [0, 0.75, 1.5]; % Scan through this to get Fig. S3d
nS = length(S); prRange = [0.1, 0.9]; nR = 100;
x = 0.01 : 0.01 : 0.99; nx = length(x);

t = 0; h = waitbar(0); lH = []; lStr = {}; figure
for iS = 1 : nS
    s = S(iS); p = 2 * s * rand(N, 1) - s + 2;
    W = randn(N) * g / sqrt(N);
    prRate = zeros(nR, 1); prCov = prRate;
    for iR = 1 : nR
        targetPR = rand*(prRange(2)-prRange(1))+prRange(1);
        [R, E] = simrepdc(W, targetMean, targetPR, p);
        prRate(iR) = pr(R); prCov(iR) = pr(E);
        t = t + 1; waitbar(t/nS/nR, h)
    end
    lH(iS) = scatter(prRate, prCov, [], co(iS, :), 'filled'); hold on
    lStr{iS} = ['width=', num2str(2*s)];

    if iS == 1
        prPred = x * (1 - g^2 * targetMean)^2;
    else
        prPred = zeros(nx, 1);
        for ix = 1 : nx
            targetPR = x(ix);
            [mu, sigma] = lognparam(targetMean, targetPR);
            fun1 = @(p) exp(2*(1-1./p)*mu+2*(1-1./p).^2*sigma^2);
            fun2 = @(p) exp(4*(1-1./p)*mu+8*(1-1./p).^2*sigma^2);
            temp1 = integral(fun1, 2-s, 2+s) / (2*s);
            temp2 = integral(fun2, 2-s, 2+s) / (2*s);
            prPred(ix) = (temp1^2/temp2) * (1 - g^2 * targetMean)^2;
        end
    end
    plot(x, prPred, 'color', co(iS, :))
end
delete(h)
xlabel('D(R)'); ylabel('D(C)'); legend(lH, lStr)

%% Decoding performance vs D(R) (Fig. 5c)
N = 100; targetMean = 10; g = 0.2;
nRepeat = 1000; % Number of repeated experiments
d0 = 1.5; % Overall scale of distance between mean responses of 2 stimuli
minRate = 0.01; % Rule out small firing rate responses
tChange = 1e4;
v0 = ones(N, 1) / sqrt(N); % Normal vector to simplex

dp = zeros(1, nRepeat); DR = dp;
h = waitbar(0);
for iRepeat = 1 : nRepeat
    % Set parameters
    targetPR = rand*0.8+0.1; DR(iRepeat) = targetPR;
    sigmaSq = -log(targetPR);
    sigma = sqrt(sigmaSq);
    mu = log(targetMean)-sigmaSq/2;
    
    % Sample distribution means to the 2 stimuli
    Z1 = lognrnd(mu, sigma, N, 1); v1 = Z1 / norm(Z1);
    Z2 = 0; t = 0;
    while min(Z2) < minRate % Must sample valid mean response vectors
        t = t + 1;
        if t > tChange % Bad Z1 sample
            Z1 = lognrnd(mu, sigma, N, 1); v1 = Z1 / norm(Z1);
            t = 0;
        end
        temp = randn(N, 1);
        % Make sure the new vector has the same mean firing rate and D(R)
        temp = temp - (v0'*temp)*v0 - (v1'*temp)*v1; 
        Z2 = Z1 + d0 * sqrt(1/targetPR-1) * temp / norm(temp);
    end

    % Generate data
    W = randn(N) * g / sqrt(N);
    temp1 = inv(diag(1./sqrt(Z1))-W); C1 = temp1 * temp1';
    temp2 = inv(diag(1./sqrt(Z2))-W); C2 = temp2 * temp2';
    
    % Estimate Fisher's linear discriminant
    d = Z2 - Z1;
    dp(iRepeat) = sqrt(2*d'*((C1+C2)\d));
    
    waitbar(iRepeat/nRepeat, h)
end
delete(h)

% Theory
x = 0.1 : 0.01 : 0.9;
dpTheory = sqrt(d0^2 * (1./x-1) .* (1./(x*targetMean)+g^2));

figure; scatter(DR, dp); hold on; plot(x, dpTheory)
xlabel('D(R)'); ylabel('d prime')
legend('Simulation', 'Theory')

