function C = generatebigc(W, Z, R)
%GENERATEBIGC Simulate a covariance matrix (before sampling and scaling).
%   Arguments:
%       W: connectivity matrix
%       Z: temporal sampling matrix
%       R: a vector of firing rates as the empirical rate distribution.
%   Returns:
%       C: the resulting large covariance matrix.
N = size(W, 1);
R = randsample(R, N, true);
Linv = 1 ./ sqrt(R);
temp = (diag(Linv) - W) \ Z;
C = temp * temp';

end