function [ npNormTimeSeries ] = npSubtract_preDF( rawTimeSeries,npTimeSeries,r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fns=fieldnames(rawTimeSeries);
cellNames=fieldnames(rawTimeSeries.(fns{1}));

for J=1:length(fns)
    for K=1:length(cellNames)
        npNormTimeSeries.(fns{J}).(cellNames{K})=rawTimeSeries.(fns{J}).(cellNames{K})-r*npTimeSeries.(fns{J}).(cellNames{K});
    end
end

end

