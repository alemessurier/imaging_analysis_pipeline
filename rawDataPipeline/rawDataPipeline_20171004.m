function [ROI_positions,rawTimeSeries,Metadata,deltaF,sampRate,fns,cellNames,Stimuli,protocol,traceByStim,cells,bsDist_meanDF,meanResponses,bootStrapResults,framesEvoked ]...
=rawDataPipeline_20171004(dir_raw,dir_processed,dir_reduced,dir_ROItemplate,StimISI,whisk,bl_length,timePostStim,timeEvoked,numBoots)

% First step in analysis of 9-whisker receptive field mapping imaging  experiments.  Can be run as a function (called by 'analysisTemplate_date') but 
% recommended to run section by section after setting directories, some variables in 'analysisTemplate_date'. AML 11/13/15 


%% register and scale raw data
cd(dir_processed)
regd=dir('*.tif');
if isempty(regd)
    ex_image=registration_res(dir_raw,dir_processed);
%   imwrite(ex_image,strcat(dir_reduced,'ex_image.tif'))
end

% %% apply a 3D median filter to registered data
% cd(dir_filt)
% filtData=dir('*.tif');
% if isempty(filtData)
%     apply_3dMed_filter(dir_processed,dir_filt)
% end

%% define ROIs 
ring=10;   gap=2;  %donut-shaped npMask beginning "gap" pixels from the outer most of ROIs and has the width of donut ring = "ring" pixels
% ring=2/gap=2 from Dan O'Connor 
% ring=10/gap=3 from Svoboda 2015
corrThreshold=0.2; %(Peron/Svoboda 2015 neuron) correlation coefficient threshold for defining correlation traces
corr=0; % 0=not consider correlation, 1=use correlation map, 2=directly check correlation with ROI signals

cd(dir_reduced);
tmp=dir('ROI_positions*');
if isempty(tmp)
    [ROI_positions,npMasks_new]=label_ROIs_HW_AML(dir_processed,[],ring,gap);
    save(strcat(dir_reduced,'ROI_positions_',date,'.mat'),'ROI_positions','npMasks_new');
else
    load([dir_reduced,tmp.name])
%     [ npMasks_new,ROI_positions] = HW_makeNpMask( ROI_positions,ring,gap);
%     save([dir_reduced,'ROI_positions.mat'],'ROI_positions','npMasks_new')
end

%% Extract raw fluorescence time series from ROIs  *****only need to run this once; saves out rawTimeSeries and Metadata to dir_reduced****

%get ROI raw time series, save
[ rawTimeSeries,Metadata ] = getFluoTimeSeries_wrapper( dir_processed,  ROI_positions);
save(strcat(dir_reduced,'Metadata'),'Metadata');
    save(strcat(dir_reduced,'ROI_timeSeries.mat'),'rawTimeSeries');

% make NPcorr masks, save    
npMasks=get_NPcorr_mask(dir_reduced,corrThreshold);
save([dir_reduced,'npMask_corr_',date,'.mat'],'npMasks');

% get npMask time series, save
[ rawTimeSeriesNP_corr,Metadata ] = getFluoTimeSeries_wrapper( dir_processed,  npMasks );
save([dir_reduced,'npCorr_TS_',date,'.mat'],'rawTimeSeriesNP_corr');
       


%% filter data with sliding average
if filtData==1
    [ filtTimeSeries ] = slidingAvg_rawF_wrapper( rawTimeSeries,ptsToAvg,type );
%     save(strcat(dir_reduced,'ROI_timeSeries.mat'),'rawTimeSeries','filtTimeSeries');
    [ npfiltTimeSeries ] = slidingAvg_rawF_wrapper( rawTimeSeriesNP_corr,ptsToAvg,type );
%     save(strcat(dir_reduced,'NP_timeSeries.mat'),'rawTimeSeriesNP','npfiltTimeSeries');
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

    [ deltaF,sampRate,truncTotal,fns,cellNames ] = calc_deltaF_wrapper(npNormTimeSeries,Metadata );


%% Create Stimuli structure containing stimulus times, based on user-entered ISI
  
 [ Stimuli,protocol ] = make_Stimuli( dir_reduced,fns,StimISI,Metadata );
% [ Stimuli,protocol,~] = HWmake_Stimuli( dir_reduced,fns,Metadata );
%% make traceByStim, find whisker tuning

framesEvoked=(ceil(bl_length*sampRate(1))+1):(ceil(bl_length*sampRate(1))+ceil(timeEvoked*sampRate(1))); 
[ traceByStim, stimFramesAll ] = make_traceByStim( whisk,Stimuli,Metadata,deltaF,bl_length,timePostStim );
[cells,~]=step2(traceByStim,framesEvoked);

%% find spontaneous activity
[ sponTrace ] = make_sponTrace( Stimuli,Metadata,deltaF,bl_length,timePostStim );

%% make labeled movie
tmp=dir(dir_processed);
num=ceil(numel(tmp)/2);
path_movie=strcat(dir_processed,tmp(num).name);
make_labeledMovie( path_movie,Stimuli,stimFramesAll,whisk,sampRate )
%% kstest for whisker responsiveness (only use if included an empty stimulus to probe spontaneous activity)

[ permTestResults ] = permuteTest_whisk_2016( sponTrace,traceByStim,numBoots,framesEvoked,permTestType );

%% save data structures

save(strcat(dir_reduced,'step1NP_',date,'.mat'),'deltaF','Stimuli',...
    'protocol','traceByStim','cells','sampRate','framesEvoked','sponTrace',...
    'permTestResults');
save([dir_reduced,'settings_',date,'mat'],'bl_length','filtData','mag',...
    'movies_exclude','numBoots','permTestType','ptsToAvg','ROIs_exclude',...
    'StimISI','timeEvoked','timePostStim','type');




