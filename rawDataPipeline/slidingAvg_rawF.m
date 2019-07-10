function [ filtTrace ] = slidingAvg_rawF( rawTrace,ptsToAvg,type )
% moving average timeseries filter. 'ptsToAvg' is number of points to
% average over; 'type' is either 'median' or 'mean', both set in
% 'analysisTemplate.m'

filtTrace=zeros(1,length(rawTrace));

for t=(1+ptsToAvg):(length(filtTrace)-ptsToAvg)
    switch type
        case 'median'
            filtTrace(t)=median(rawTrace((t-ptsToAvg):(t+ptsToAvg)));
        case 'mean'
            filtTrace(t)=mean(rawTrace((t-ptsToAvg):(t+ptsToAvg)));
    end
end
end

