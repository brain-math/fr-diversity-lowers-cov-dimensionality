function E = repeatedsample(W, Z, R, nRepeat)
%REPEATEDSAMPLE Aggregate the eigenvalues of nRepeat submatrices.
N = size(W, 1); S = length(R);
C = generatebigc(W, Z, R);
E = zeros(nRepeat*S, 1);
for iRepeat = 1 : nRepeat
    observed = randperm(N, S);
    Ch = C(observed, observed);
    try
        E((iRepeat-1)*S+1:iRepeat*S) = svd(Ch);
    catch
        E((iRepeat-1)*S+1:iRepeat*S) = 0;
    end
end

end