% Numerical verification of the main result of Tian et al. (2024)
% Code by: Gengshuo John Tian

set(groot, 'defaultLineLineWidth', 2)
set(groot, 'defaultAxesLineWidth', 2)
set(groot, 'defaultAxesFontSize', 20)
set(groot, 'defaultAxesTickDir', 'out')
set(groot,  'defaultAxesTickDirMode', 'manual')

%% Fig. 2a (left)
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
