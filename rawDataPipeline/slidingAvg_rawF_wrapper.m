function [ filtTimeSeries ] = slidingAvg_rawF_wrapper( rawTimeSeries,ptsToAvg,type )
%wrapper fxn for moving average timeseries filter. 'ptsToAvg' is number of points to
% average over; 'type' is either 'median' or 'mean', both set in
% 'analysisTemplate.m'.

fns=fieldnames(rawTimeSeries);
cellNames=fieldnames(rawTimeSeries.(fns{1}));

for K=1:length(fns)
    for J=1:length(cellNames)
        filtTimeSeries.(fns{K}).(cellNames{J})=slidingAvg_rawF(rawTimeSeries.(fns{K}).(cellNames{J}),ptsToAvg,type);
    end
end

end

