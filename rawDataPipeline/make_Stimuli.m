function [ Stimuli,protocol ] = make_Stimuli( dir_reduced,fns,StimISI,Metadata )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

tmp=dir(fullfile(dir_reduced,'/*.ibt'));
stimFiles=arrayfun(@(x)x.name,tmp,'Uni',0);
stimFileNums=cellfun(@(x)x((strfind(x,'.')-2):(strfind(x,'.')-1)),stimFiles,'Uni',0);
for i=1:length(stimFileNums)
    repl=strfind(stimFileNums{i},'_');
    if ~isempty(repl)
        stimFileNums{i}(repl)='0';
    else
    end
end

for K=1:length(fns)
    fn=fns{K};
    fn_num=fn(end-1:end);
         sampRate=1/(Metadata.(fn).acqNumAveragedFrames*Metadata.(fn).acqScanFramePeriod);


    if ismember(fn_num,stimFileNums)
        protocol.(fn)='stim';
        clear data
        clear param
        clear stimTimes
        
        thisFileInds=strcmp(stimFileNums,fn_num);

        ibt_path=strcat(dir_reduced,stimFiles{thisFileInds});
        [data,param]=ibtRead(ibt_path);
        realISIs{K}=diff(param.sweepTimes);
        whiskOrder{K}=param.stim.piezoNumber;
        
        % find stim times for aligning to imaging data
        
        sweepNum=1;
        PrestimDelay=param.stim.preDelay(sweepNum)/1000;
   
        
        nstimuli = floor((Metadata.(fn).NumFrames/sampRate) / StimISI);    % enough stimuli to fill the movie
    Stimuli.(fn).IgorSweepNo = 0:1:(nstimuli-1);      % Igor sweep number starts at zero.
    Stimuli.(fn).Label = param.stim.num(1:nstimuli);      
    Stimuli.(fn).Time = zeros([nstimuli 1]);
    Stimuli.(fn).Time = Stimuli.(fn).IgorSweepNo * StimISI + PrestimDelay;   % time in sec from initial trigger to each stimulus
    
    % Register Stimulus times to the movie, with no initial blackout period
%     [Stimuli.(fn), Metadata.(fn)] = RegisterStimuliToMovie(Stimuli.(fn), Metadata.(fn), 0)
    else
        protocol.(fn)='spon';
    end

end


end

