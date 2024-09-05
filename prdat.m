function d = prdat(eigenvalues, alpha, nOl)
%PRDAT Calculate participation ratio from the eigenvalue distribution.
% Corrected for temporal sampling and outliers.
% alpha: N / T
if max(eigenvalues) <= 0
    d = 0;
    disp('Illegitimate eigenvalues!')
    return
end
eigenvalues = sort(eigenvalues, 'descend');
eigenvalues = eigenvalues(nOl+1:end);
d = pr(eigenvalues);
d = dhattod(d, alpha); % Account for temporal sampling

end

