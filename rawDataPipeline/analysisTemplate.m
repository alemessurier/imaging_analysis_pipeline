% analysis template - calls 'reduce_image_data'  
% Experiment: ___ (mouse), imaged on ___
% 
% this file should be renamed for each experiment and saved in dir_reduced 
% (set below), along with all igor binary files (.ibt) with stimulus info
%% set directories, define some parameters

dir_raw='E:\Data\raw\DENR9a\20161021\f1\to_reg\';   % set directory containing raw tiff files collected w/scanImage
dir_processed='E:\Data\raw\DENR9a\20161021\f1\registered\';    % set directory to save registered data 
dir_reduced='E:\Data\reduced\DENR9a\20161021\f1\';  % set directory for reduced data
dir_ROItemplate=[]; % optional: if same field was imaged previously, directory to look for image masks


bl_length=0.5; % pre-stim baseline in seconds
timePostStim=3; % time to include post-stim in 'traceByStim'; seconds
timeEvoked=1; % time post-stimulus to include for measures of whisker response; seconds 
StimISI=5.01; % seconds
whisk={'e3' 'e2' 'e1' 'd3' 'd2' 'd1' 'c3' 'c2' 'c1' }; % set of whiskers deflected, by row starting from most dorsal and caudal whisker.
numBoots=10000; % number of iterations for permutation tests of whisker responsiveness
filtData=1; % if 1, ROI raw fluorescence timeseries will be filtered with a moving average filter; if 0, no temporal filtering
ptsToAvg=2; % number of frames before and after each frame to average over for temporal filter
filter_type='median'; % type of temporal filter (can be 'median' or 'mean')
depth=140; % depth of imaging field from dura
mag=1.5; % software zoom set in ScanImage during imaging
permTestType='median'; % type of permutation test to use for testing for whisker responsiveness
ring=10;   gap=2;  %donut-shaped npMask beginning "gap" pixels from the outer most of ROIs and has the width of donut ring = "ring" pixels
% ring=2/gap=2 from Dan O'Connor 
% ring=10/gap=3 from Svoboda 2015
corrThreshold=0.2; %(Peron/Svoboda 2015 neuron) correlation coefficient threshold for defining correlation traces

%% set matlab path

addpath(genpath('E:\imaging_analysis_pipeline\rawDataPipeline\')); % directory containing functions for analysis pipeline
addpath(genpath('E:\imaging_analysis_pipeline\plotting_2017\')); % directory containing plotting tools
addpath(dir_raw,dir_processed,dir_reduced); % add data directories
if ~isempty(dir_ROItemplate)
    addpath(dir_ROItemplate);
end

%% run analysis pipeline to reduce and save data

rawDataPipeline_shared(dir_raw,dir_processed,dir_reduced,filtData,...
    filter_type,ptsToAvg,mag,StimISI,whisk,bl_length,timePostStim,timeEvoked,numBoots,permTestType)