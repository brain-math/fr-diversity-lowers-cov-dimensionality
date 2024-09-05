function [E, R, T] = gatherer(areaData)
%GATHERER Gather all E and R for one areaData into cell arrays
nChunk = length(areaData);
for iChunk = 1 : nChunk
    E{iChunk} = areaData(iChunk).E;
    R{iChunk} = areaData(iChunk).R;
    T{iChunk} = size(areaData(iChunk).X, 1);
end

end