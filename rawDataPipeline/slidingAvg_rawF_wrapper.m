function [ filtTimeSeries ] = slidingAvg_rawF_wrapper( rawTimeSeries,ptsToAvg,type )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

fns=fieldnames(rawTimeSeries);
cellNames=fieldnames(rawTimeSeries.(fns{1}));

for K=1:length(fns)
    for J=1:length(cellNames)
        filtTimeSeries.(fns{K}).(cellNames{J})=slidingAvg_rawF(rawTimeSeries.(fns{K}).(cellNames{J}),ptsToAvg,type);
    end
end

end

