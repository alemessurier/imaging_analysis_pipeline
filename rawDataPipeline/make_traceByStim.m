function [ traceByStim,stimFramesAll ] = make_traceByStim( whisk,Stimuli,Metadata,deltaF,bl_length,timePostStim )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
fns=fieldnames(Stimuli);
cellNames=fieldnames(deltaF.(fns{1}));

for i=1:length(cellNames)
    cn=cellNames{i};
    for j=1:length(whisk)
        whisker=whisk{j};
        traceByStim.(cn).(whisker)=[];
    end
end


for K=1:length(fns)
    fn=fns{K};
     sampRate=1/(Metadata.(fn).acqNumAveragedFrames*Metadata.(fn).acqScanFramePeriod);
    truncTotal = Metadata.(fn).NumFrames-length(deltaF.(fn).(cellNames{1}));
    
%      if strcmp(protocol.(fn),'stim')==1;
            
        % find stim times for aligning to imaging data
          stimFrames=floor(Stimuli.(fns{K}).Time*sampRate);%-truncTotal);
          stimOrder=Stimuli.(fn).Label(stimFrames(1:(end-1))>0);
        stimFrames=stimFrames(stimFrames>0);
       stimFramesAll.(fn)=stimFrames;
        
        % index stimtimes by whisker identity
        
        bl_im=ceil(bl_length*sampRate); %pre-stim baseline in frames
        frames_postStim=ceil(timePostStim*sampRate); %post-stim period to include
         
        for i=1:length(whisk);
            whiskNum=i-1;
            stimInds{i}=stimOrder==whiskNum;
        end
        
        stimFrames_whisk=cellfun(@(x)stimFrames(x),stimInds,'Uni',0);
        
        
        for k=1:length(whisk)
            stimFrames_thisWhisk=stimFrames_whisk{k};
            stimFrames_thisWhisk=stimFrames_thisWhisk(stimFrames_thisWhisk>ceil(bl_im) & stimFrames_thisWhisk<(length(deltaF.(fn).(cellNames{1}))-frames_postStim));% bug fixed 5/9/16
            whisker=whisk{k};
            for j=1:length(cellNames)
                cn=cellNames{j};
                %baseline subtraction:
                stimBlock=arrayfun(@(x)deltaF.(fn).(cn)((x-bl_im):(x+frames_postStim))-mean(deltaF.(fn).(cn)((x-(bl_im)):x)),stimFrames_thisWhisk,'Uni',0); 
                traceByStim.(cn).(whisker)=[traceByStim.(cn).(whisker); stimBlock'];
            end
        end
%     else
        
%     end
    
end


for j=1:length(cellNames)
    cn=cellNames{j};
    for i=1:length(whisk)
        whisker=whisk{i};
        traceByStim.(cn).(whisker)=vertcat(traceByStim.(cn).(whisker){:});
    end
end



end

