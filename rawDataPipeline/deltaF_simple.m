function [ deltaF ] = deltaF_simple( trace)
% New version 8/5/15, post-thesis commitee meeting

percentile_baseline=0.05;        
        bltmp=sort(trace,'ascend');
        bline=bltmp(floor(percentile_baseline*length(bltmp)));
        if le(bline,0)
            findBl=bltmp(bltmp>0);
            if ~isempty(findBl)
                
                bline=findBl(1);
            end
        end
        
        
        
   
deltaF=(trace - bline)/bline;









end

