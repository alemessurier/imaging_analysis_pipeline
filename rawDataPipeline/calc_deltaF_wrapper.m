function [ deltaF,sampRate,fns,cellNames ] = calc_deltaF_wrapper(rawTimeSeries,Metadata )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fns=fieldnames(rawTimeSeries);
cellNames=fieldnames(rawTimeSeries.(fns{1}));

    
for J=1:length(fns)
    fn=fns{J};
    
    sampRate(J)=1/(Metadata.(fn).acqNumAveragedFrames*Metadata.(fn).acqScanFramePeriod);
    orig_recorded_frames(J) = Metadata.(fn).NumFrames;
    recordedSeconds(J)=orig_recorded_frames(J)/sampRate(J);

    for K=1:length(cellNames)
        cn=cellNames{K};
        [ deltaF.(fn).(cn)] = svoboda_deltaF_simple( rawTimeSeries.(fn).(cn),0 );

    end
end



end

