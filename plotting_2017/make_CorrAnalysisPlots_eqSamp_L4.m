function  [fullfield_NH,fullfield_EN,inBarr_NH,inBarr_EN,xBarr_NH,xBarr_EN]=make_CorrAnalysisPlots_eqSamp_L4( paths_NH,paths_EN,r,ROIs_NH,ROIs_EN )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% only use fields that have at least 3 ROIs 
useEN=cellfun(@(x)length(x)>2,ROIs_EN);
paths_EN=paths_EN(useEN);
ROIs_EN=ROIs_EN(useEN);

useNH=cellfun(@(x)length(x)>2,ROIs_NH);
paths_NH=paths_NH(useNH);
ROIs_NH=ROIs_NH(useNH);

%% full field



for K=1:length(paths_EN)
             [~,traceByStim,sponTrace,framesEvoked,~,...
                ~,ROIsInBarrel,~,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L4( paths_EN{K},r );
                 traceByStim=traceByStim(2);
                    sponTrace=sponTrace(2);
           
       
   
    sigROIs=ROIs_EN{K};
    
    results_field_EN(K) = find_corr_simple( traceByStim,sponTrace,sigROIs,framesEvoked,ROI_positions,mag);
    
end

%% in field

for K=1:length(paths_NH)
    
    
            [~,traceByStim,sponTrace,framesEvoked,~,...
                ~,ROIsInBarrel,~,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L4( paths_NH{K},r );
                 
    sigROIs=ROIs_NH{K};
    results_field_NH(K) = find_corr_simple( traceByStim,sponTrace,sigROIs,framesEvoked,ROI_positions,mag);
    
end
%% concatenate results across fields

NC_NH_all=cat(1,results_field_NH(:).NC);
NC_NH_all=mean(NC_NH_all,2);

ROIdist_NH_all=cat(2,results_field_NH(:).ROIdistance);
SC_NH_all=cat(1,results_field_NH(:).SC);

fullfield_NH.NC=NC_NH_all;
fullfield_NH.SC=SC_NH_all;
fullfield_NH.dist=ROIdist_NH_all;


NC_EN_all=cat(1,results_field_EN(:).NC);
NC_EN_all=mean(NC_EN_all,2);

ROIdist_EN_all=cat(2,results_field_EN(:).ROIdistance);
SC_EN_all=cat(1,results_field_EN(:).SC);

fullfield_EN.NC=NC_EN_all;
fullfield_EN.SC=SC_EN_all;
fullfield_EN.dist=ROIdist_EN_all;


%% make plots comparing SC,NC by distance within column pairs vs all

% find NC,SC for ROIs within and across columns
results_NH=[];
for K=1:length(paths_NH)
   
        
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L4( paths_NH{K},r );
           
    sigROIs=ROIs_NH{K};
     [ results ] = find_corrInBarrelvsAcross( sigROIs,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
     if ~isempty(results)
        results_NH=[results_NH,results];
     end
end

results_EN=[];
for K=1:length(paths_EN)
 
        
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L4( paths_EN{K},r );
            
    sigROIs=ROIs_EN{K};
    [ results ] = find_corrInBarrelvsAcross( sigROIs,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
    if ~isempty(results)
        results_EN=[results_EN,results];
    end
end





%%
  NC_EN_inBarr=cat(1,results_EN(:).NC_inBarr);
inBarr_EN.NC=mean(NC_EN_inBarr,2);
NC_EN_xBarr=cat(1,results_EN(:).NC_xBarr);
xBarr_EN.NC=mean(NC_EN_xBarr,2);
inBarr_EN.SC=cat(1,results_EN(:).SC_inBarr);
xBarr_EN.SC=cat(1,results_EN(:).SC_xBarr);
inBarr_EN.dist=cat(2,results_EN(:).ROIdist_inBarr);
xBarr_EN.dist=cat(2,results_EN(:).ROIdist_xBarr);



 NC_NH_inBarr=cat(1,results_NH(:).NC_inBarr);
inBarr_NH.NC=mean(NC_NH_inBarr,2);
NC_NH_xBarr=cat(1,results_NH(:).NC_xBarr);
xBarr_NH.NC=mean(NC_NH_xBarr,2);
inBarr_NH.SC=cat(1,results_NH(:).SC_inBarr);
xBarr_NH.SC=cat(1,results_NH(:).SC_xBarr);
inBarr_NH.dist=cat(2,results_NH(:).ROIdist_inBarr);
xBarr_NH.dist=cat(2,results_NH(:).ROIdist_xBarr);



end