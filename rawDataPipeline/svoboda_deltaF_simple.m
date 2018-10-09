function [ deltaF] = svoboda_deltaF_simple( trace,kurt_thresh )
% New version 8/5/15, post-thesis commitee meeting

kurt=kurtosis(trace);
percentile_baseline=0.05;

r=[];
samples=length(trace);

switch kurt>kurt_thresh
    case 1
        
        bltmp=sort(trace,'ascend');
        bline=bltmp(floor(percentile_baseline*length(bltmp)));
        if le(bline,0)
            findBl=bltmp(bltmp>0);
            if ~isempty(findBl)
                
                bline=findBl(1);
            end
        end
        
        
        
    case 0
        bline=median(trace);
end

for t = 1:samples
    r = [r; (trace(t) - bline)/bline];
end






deltaF=r;

end

