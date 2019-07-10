function [ sponTrace ] = make_sponTrace( Stimuli,sampRate,deltaF,bl_length,timePostStim )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

fns=fieldnames(Stimuli);
cellNames=fieldnames(deltaF.(fns{1}));

for i=1:length(cellNames)
    cn=cellNames{i};
    sponTrace.(cn)=[];
end


for K=1:length(fns)
    fn=fns{K};
%     sampRate=1/(Metadata.(fn).acqNumAveragedFrames*Metadata.(fn).acqScanFramePeriod);
%     truncTotal = Metadata.(fn).NumFrames-length(deltaF.(fn).(cellNames{1}));
%     
%     
%     
     % find stim times for aligning to imaging data
    stimFrames=floor(Stimuli.(fns{K}).Time*sampRate);%-truncTotal);
    stimOrder=Stimuli.(fn).Label(stimFrames(1:(end-1))>0);
    stimFrames=stimFrames(stimFrames>0);
        
    % index stimtimes by whisker identity
    
    bl_im=ceil(bl_length*sampRate); %pre-stim baseline in frames
    frames_postStim=ceil(timePostStim*sampRate); %post-stim period to include
    
    sponInds=stimOrder==9; %change back to 9!!
    
    sponFrames=stimFrames(sponInds);
    
    
    sponFrames=sponFrames(sponFrames>ceil(bl_im) & sponFrames<(length(deltaF.(fn).(cellNames{1}))-frames_postStim));
    
    for j=1:length(cellNames)
        cn=cellNames{j};
        %baseline subtraction:
        sponBlock=arrayfun(@(x)deltaF.(fn).(cn)((x-bl_im):(x+frames_postStim))-mean(deltaF.(fn).(cn)((x-(bl_im)):x)),sponFrames,'Uni',0);
%         sponBlock=arrayfun(@(x)deltaF.(fn).(cn)((x+bl_im):(x+frames_postStim))-mean(deltaF.(fn).(cn)(x:(x+bl_im))),sponFrames,'Uni',0); %2/2/18 shifted window later in trace

        sponTrace.(cn)=[sponTrace.(cn); sponBlock'];
    end
end




for j=1:length(cellNames)
    cn=cellNames{j};
    sponTrace.(cn)=vertcat(sponTrace.(cn){:});
end



end

