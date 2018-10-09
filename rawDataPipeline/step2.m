function [cells,Wh ] = step2( traceByStim,inds )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cellNames=fieldnames(traceByStim);
whisk=fieldnames(traceByStim.(cellNames{1}));

for K=1:numel(cellNames)
    cn=cellNames{K};
%     thresh=stats.(cn).meanDFall+(2*stats.(cn).stdevDFall);
    for W=1:numel(whisk)
        whisker=whisk{W};
        intdF=trapz(traceByStim.(cn).(whisker)(:,(inds)),2);%find area under curve for each stimulus presentation
        meanIntdF(W)=mean(intdF);
        medianIntdF(W)=median(intdF);
        %calculate top 30% of responses
        intdF=sort(intdF,'descend');
        top30=intdF(1:floor(0.3*length(intdF)));
        meanTop30(W)=mean(top30);
        
        %calculate values of max +/- 1 frames
        dFs=mean(traceByStim.(cn).(whisker)(:,inds),2);
        maxes=max(traceByStim.(cn).(whisker)(:,inds),[],2);
        meanPeakDF(W)=mean(maxes);
        meanDF(W)=mean(dFs);
        medianResp=median(traceByStim.(cn).(whisker),1);
        medianDF(W)=mean(medianResp(inds));
        
     
        
        meanRWhiskbyCell(W,K)=meanIntdF(W);
        Wh.(whisker).meanIntdF(K)=meanIntdF(W); %make structure of mean responses to this whisker for all cells
        %determine P(r)for each whisker using threshold crossings
%         stimEpochs=traceByStim.(cn).(whisker)(:,inds);
%         tmp=stimEpochs>=thresh;
%         tmp=sum(tmp,2);       
%         Pr(W)=(sum(tmp~=0))/(length(tmp));
%         Wh.(whisker).Pr(K)=Pr(W);
%         PrWhiskByCell(W,K)=Pr(W);
%         cells.Prs.(cn).(whisker)=Pr(W);
        cells.AUCs.(cn).(whisker)=meanIntdF(W);
    end
    cells.(cn).meanTop30=meanTop30;
    cells.(cn).meanIntdFWhisk=meanIntdF; %for each cell, save mean response to each whisker in a vector
     cells.(cn).medianIntdFWhisk=medianIntdF;
    [PWr(K),PWinds]=max(meanIntdF); %find PW for each cell and response size to PW
    cells.(cn).meanDF=meanDF;
    cells.(cn).meanPeakDF=meanPeakDF;
    cells.(cn).medianDF=medianDF;
    PW_response{K}=whisk{PWinds}; 
    cells.(cn).PW_response=PW_response{K}; %save PW and PWr for this cell
    cells.(cn).PW_response_mag=PWr(K);
    meanDiffPWr(K)=sum(abs(meanIntdF-PWr(K)))/8; %compute difference between response to PW and response to other whiskers
    cells.(cn).meanDiffPWr=meanDiffPWr(K);
%     cells.(cn).Pr_whisk=Pr;
    cells.(cn).whisk=whisk;
%     [PrToPW(K),PW_Pr_inds]=max(Pr);
%     PW_Pr{K}=whisk{PW_Pr_inds};
%     cells.(cn).PW_Pr=PW_Pr{K};
%     cells.(cn).PrToPW=PrToPW(K);
end

%Sort cells in descending order according to PrToPW
%[~,inds]=sort(PrToPW,'descend');
[~,inds]=sort(PWr,'descend');
cells.cellsSorted=cellNames(inds);


