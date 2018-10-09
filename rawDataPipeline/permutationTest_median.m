function [ P ] = permutationTest_median(stimLocked,spont,numReps )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

T=abs(nanmedian(stimLocked)-nanmedian(spont));
dim1=find(size(stimLocked)>1);
pooled=cat(dim1,stimLocked,spont);

for t=1:numReps
    indsA=randsample(1:length(pooled),length(stimLocked));
    A=pooled(indsA);
    indsB=ones(length(pooled),1);
    indsB(indsA)=0;
    indsB=logical(indsB);
    B=pooled(indsB);    
    PT(t)=abs(nanmedian(A)-nanmedian(B));
end
P=sum(PT>T)/length(PT);



end

