function [DR, DC] = getpr(areaData, nOl)
%GETPR Calculate D(C) and D(R) after eliminating outliers
nChunk = length(areaData); DR = zeros(1, nChunk); DC = zeros(1, nChunk);
for iChunk = 1 : nChunk
    [T, N] = size(areaData(iChunk).X); alpha = N / T;
    DR(iChunk) = pr(areaData(iChunk).R);
    DC(iChunk) = prdat(areaData(iChunk).E, alpha, nOl);
end

end