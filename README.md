# imaging_analysis_pipeline
MATLAB functions for processing 2photon calcium imaging data into reduced data for analysis of 9-whisker receptive fields.

## basics
Main script 'rawDataPipeline_shared' calls functions for reducing raw calcium imaging movie files collected with ScanImage. First, use 'analysisTemplate' to manually enter some variables, including directories containing raw data and stimulus files, directories for saving reduced data, and some parameters for data analysis. 'analysisTemplate' must be saved in 'dir_reduced'. Running 'analysisTemplate' will call 'rawDataPipeline' as a function; however, I recommend running the first two sections of 'analysisTemplate', then running 'rawDataPipeline' section-by-section since some image processing steps may take a long time, and data from each step is saved incrementally.

## data
Raw movie files must be 512x512 pixel tiff stacks. These must be saved in directory 'dir_raw' defined in 'analysisTemplate'. Registered movies will be saved in 'dir_processed', and reduced data will be saved to 'dir_reduced'. Igor binary files containing stimulus information (.ibt) for each movie must be saved in 'dir_reduced', and names must match movie names.

## reduced data outputs

'ROI_positions_date.mat':
 - ROI_positions: 512x512xn logical array of ROI masks; n is number of ROIs circled
 - npMasks_basic: 512x512xn logical array of neuropil masks for each ROI
 
 'Metadata.mat': metadata for each movie extracted from tiff headers
 
 'ROI_timeSeries.mat':
  - rawTimeSeries: structure with fields for each ROI in each movie; raw fluorescence time series by frame
  - filtTimeSeries: moving-average filtered fluorescence time series

'npCorr_TS_date.mat': 
 - rawTimeSeriesNP_corr: time series for neuropil masks
 - filtTimeSeries: moving-average filtered fluorescence time series
 
'step1_date.mat':
 - deltaF: structure of deltaF/F time series for each ROI for each movie.
 - Stimuli: structure of stimulus information for each movie, including whisker ID and frames aligned to stim delivery.
 - traceByStim: structure of dF/F in time window before and after each stim iteration by ROI and by whisker identity.
 - sampRate: vector of frame rates for each movie
 - framesEvoked: vector of frames after each stimulus to measure evoked dF/F
 - sponTrace: structure of dF/F for each ROI on blank trials
 - permTestResults: permutation test derived p-values for each ROI comparing whisker evoked dF/F to blank trial dF/F.
 
 'settings_date.mat': contains variables set in 'analysisTemplate.m'
