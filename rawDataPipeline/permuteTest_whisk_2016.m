function [ pvals ] = permuteTest_whisk_2016( sponTrace,traceByStim,numReps,framesEvoked,type )
%for each ROI, tests for difference in means between spontaeous dF/F and
%dF/F evoked by each whisker. Returns p-value for each whisker. For further analysis must be corrected for 9 comparisons using FDR. 
%'numReps' and 'type' are set in 'analysisTemplate.m'

cellNames=fieldnames(sponTrace);
whisk=fieldnames(traceByStim.(cellNames{1}));

for K=1:length(cellNames)
    sponMeans=mean(sponTrace.(cellNames{K})(:,framesEvoked),2);
    cn=cellNames{K};
    
    for J=1:length(whisk)
        whiskMeans{J}=mean(traceByStim.(cn).(whisk{J})(:,framesEvoked),2);
    end
    
    switch type
        case 'mean'
            
            parfor J=1:length(whisk)
                tmpPVal(J)=permutationTest(whiskMeans{J},sponMeans,numReps );
            end
            
        case 'median'
            
            parfor J=1:length(whisk)
                tmpPVal(J)=permutationTest_median(whiskMeans{J},sponMeans,numReps );
            end
    end
    
    for J=1:length(whisk)
        pvals.(cn).(whisk{J}) = tmpPVal(J);
    end
    
end
end

