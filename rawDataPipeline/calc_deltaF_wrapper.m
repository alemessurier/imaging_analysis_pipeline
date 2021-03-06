function [ deltaF,sampRate,truncTotal,fns,cellNames ] = calc_deltaF_wrapper(rawTimeSeries,Metadata )
%Calculates deltaF/F for each ROI in each movie
fns=fieldnames(rawTimeSeries);
cellNames=fieldnames(rawTimeSeries.(fns{1}));
kurt=zeros(length(fns),length(cellNames));
    
for J=1:length(fns)
    fn=fns{J};
    
    sampRate(J)=1/(Metadata.(fn).acqNumAveragedFrames*Metadata.(fn).acqScanFramePeriod);
    orig_recorded_frames(J) = Metadata.(fn).NumFrames;
    recordedSeconds(J)=orig_recorded_frames(J)/sampRate(J);

    for K=1:length(cellNames)
        cn=cellNames{K};
         [ deltaF.(fn).(cn)] = deltaF_simple( rawTimeSeries.(fn).(cn));

    end
     truncTotal(J) = orig_recorded_frames(J)-length(deltaF.(fn).(cellNames{1}));
end



end

