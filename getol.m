function nOl = getol(areaData, fitted, fitParam)
%GETOL Determine the number of outliers in an area
nChunk = length(areaData);
nOlAll = zeros(1, nChunk);
for iChunk = 1 : nChunk
    nOlAll(iChunk) = detectol(areaData(iChunk).E, areaData(iChunk).R,...
        size(areaData(iChunk).X, 1), fitted, fitParam.nRepeat, fitParam.p);
end
% Assume a consistent number of outliers is present for each condition
nOl = max(nOlAll);

end