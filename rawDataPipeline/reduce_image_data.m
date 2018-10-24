function reduce_image_data(dir_raw,dir_processed,dir_reduced,dir_ROItemplate,...
    StimISI,whisk,bl_length,timePostStim,timeEvoked,numBoots,ring,gap,corrThreshold)

% First step in analysis of 9-whisker receptive field mapping imaging  experiments.  Can be run as a function (called by 'analysisTemplate_date') but 
% recommended to run section by section after setting directories, some variables in 'analysisTemplate_date'. AML 11/13/15 

%% register and scale raw data

cd(dir_processed)
regd=dir('*.tif');
if isempty(regd)
   registration_res(dir_raw,dir_processed);
end


%% define ROIs 

% check to see if ROI masks have already been set
cd(dir_reduced);
exist_ROIpos=dir('ROI_positions*');

if isempty(exist_ROIpos) % if ROI_positions can't be found, run label_ROIs GUI
    [ROI_positions,npMasks]=label_ROIs(dir_processed,[],ring,gap,5);
    save([dir_reduced,'ROI_positions_',date,'.mat'],'ROI_positions','npMasks');
else
    try % try loading ROI and NP masks
        load([dir_reduced,exist_ROIpos.name],'ROI_positions','npMasks')
    catch % if npMasks doesn't exist, just run function for making npMasks
        warning('npMasks does not exist, running HW_makeNpMask')
        [ npMasks,ROI_positions] = HW_makeNpMask( ROI_positions,ring,gap);
        save([dir_reduced,'ROI_positions_',date,'.mat'],'ROI_positions','npMasks');
    end
end

%% Extract raw fluorescence time series from ROIs  *****only need to run this once; saves out rawTimeSeries and Metadata to dir_reduced****

try 
    load([dir_reduced,'ROI_timeSeries.mat'],'rawTimeSeries')
    load([dir_reduced,'Metadata.mat'],'Metadata')
catch ME
    warning('rawTimeSeries does not exist - running getFluoTimeSeries')
    
    [ rawTimeSeries,Metadata ] = getFluoTimeSeries_wrapper( dir_processed,  ROI_positions);
    save(strcat(dir_reduced,'Metadata.mat'),'Metadata');
    save(strcat(dir_reduced,'ROI_timeSeries.mat'),'rawTimeSeries');
end
    
%get ROI raw time series, save
[ rawTimeSeries,Metadata ] = getFluoTimeSeries_wrapper( dir_processed,  ROI_positions);
save(strcat(dir_reduced,'Metadata'),'Metadata');
    save(strcat(dir_reduced,'ROI_timeSeries.mat'),'rawTimeSeries');

% make NPcorr masks, save    
try 
    load([dir_reduced,'NP_timeSeries.mat'],'rawTimeSeriesNP')
catch ME
    cd(dir_reduced);
    exist_npMaskCorr=dir('npMask_corr*');
    if ~isempty(exist_npMaskCorr)
        load([dir_reduced,exist_npMaskCorr.name])
    else
        % make NPcorr masks, save    
        npMasks=get_NPcorr_mask(dir_reduced,corrThreshold);
        save([dir_reduced,'npMask_corr_',date,'.mat'],'npMasks');
    end
    
    rawTimeSeriesNP = getFluoTimeSeries_wrapper( dir_processed,  npMasks );
    save([dir_reduced,'NP_timeSeries.mat'],'rawTimeSeriesNP');
end
       


%% filter data with sliding average
if filtData==1
    [ filtTimeSeries ] = slidingAvg_rawF_wrapper( rawTimeSeries,ptsToAvg,type );
    [ npfiltTimeSeries ] = slidingAvg_rawF_wrapper( rawTimeSeriesNP,ptsToAvg,type );
end

%% look through raw fluorescence movies to remove movies with large change in fluorescence

[filtTimeSeries,movies_exclude]=check_rawFluor( filtTimeSeries );
 npfiltTimeSeries=rmfield(npfiltTimeSeries,movies_exclude);
% filtTimeSeries=rmfield(filtTimeSeries,movies_exclude);
% 


%% remove ROIs with overcorrection from np subtraction
[filtTimeSeries,npfiltTimeSeries,ROIs_exclude]=check_rawROIs( filtTimeSeries,npfiltTimeSeries );
%% pre deltaF neuropil subtraction

[ npNormTimeSeries ] = npSubtract_preDF( filtTimeSeries,npfiltTimeSeries,0.3 );
%% calculate deltaF/F

    [ deltaF,sampRate,fns,cellNames ] = calc_deltaF_wrapper(npNormTimeSeries,Metadata );


%% Create Stimuli structure containing stimulus times, based on user-entered ISI
  
 Stimuli = make_Stimuli( dir_reduced,fns,StimISI,Metadata );

%% make traceByStim, find whisker tuning

framesEvoked=(ceil(bl_length*sampRate(1))+1):(ceil(bl_length*sampRate(1))+ceil(timeEvoked*sampRate(1))); 
traceByStim = make_traceByStim( whisk,Stimuli,Metadata,deltaF,bl_length,timePostStim );


%% find spontaneous activity
[ sponTrace ] = make_sponTrace( Stimuli,Metadata,deltaF,bl_length,timePostStim );

%% kstest for whisker responsiveness (only use if included an empty stimulus to probe spontaneous activity)

[ permTestResults ] = permuteTest_whisk_2016( sponTrace,traceByStim,numBoots,framesEvoked,permTestType );

%% save data structures

save(strcat(dir_reduced,'step1NP_',date,'.mat'),'deltaF','Stimuli',...
    'protocol','traceByStim','sampRate','framesEvoked','sponTrace',...
    'permTestResults');
save([dir_reduced,'settings_',date,'mat'],'bl_length','filtData','mag',...
    'movies_exclude','numBoots','permTestType','ptsToAvg','ROIs_exclude',...
    'StimISI','timeEvoked','timePostStim','type');




