function  [fullfield_NH,fullfield_EN,inBarr_NH,inBarr_EN,xBarr_NH,xBarr_EN,SW_NH,SW_EN]=make_CorrAnalysisPlots_eqSamp( paths_NH,paths_EN,npSub,ROIs_NH,ROIs_EN )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% only use fields that have at least 3 ROIs
if ~isempty(ROIs_EN)
    useEN=cellfun(@(x)length(x)>2,ROIs_EN);
    paths_EN=paths_EN(useEN);
    ROIs_EN=ROIs_EN(useEN);
end

if ~isempty(ROIs_NH)
    useNH=cellfun(@(x)length(x)>2,ROIs_NH);
    paths_NH=paths_NH(useNH);
    ROIs_NH=ROIs_NH(useNH);
end

%% full field
f1=figure; hold on
f2=figure; hold on
cmap=brewermap(length(paths_EN),'Set1');

for K=1:length(paths_EN)
    switch npSub
        case 0
            [~,traceByStim,sponTrace,framesEvoked,permTestResults,...
                ~,ROIsInBarrel,~,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_EN{K} );
        case 1
            [~,traceByStim,sponTrace,framesEvoked,permTestResults,...
                ~,ROIsInBarrel,~,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_EN{K} );
            traceByStim=traceByStim(2);
            sponTrace=sponTrace(2);
            permTestResults=permTestResults(2);
            
    end
    
    if ~isempty(ROIs_EN)
        sigROIs=ROIs_EN{K};
    else
        sigROIs=find_sigROIs(permTestResults,traceByStim);
    end
    
    results_field_EN(K) = find_corr_simple( traceByStim,sponTrace,sigROIs,framesEvoked,ROI_positions,mag);
%     [ results_SW_EN(K) ] = find_corrToSWs(sigROIs,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
%     
%     if ~isempty(results_SW_EN(K).ROIdist_inArc)
%          NC_plot=mean(results_SW_EN(K).NC_inArc,2);
%         figure(f1)
%         plot_scatterRLine_color(results_SW_EN(K).ROIdist_inArc,NC_plot,cmap(K,:))
%     end
    
  
end


figure(f1)
xlabel('distance, um')
        ylabel('NC (in arc, EN)')

        figure(f2)
xlabel('distance, um')
        ylabel('NC (in arc, NH)')

        %% in field
cmap=brewermap(length(paths_NH),'Set1');
for K=1:length(paths_NH)
    switch npSub
        case 0
            [~,traceByStim,sponTrace,framesEvoked,permTestResults,...
                ~,ROIsInBarrel,~,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_NH{K} );
        case 1
            [~,traceByStim,sponTrace,framesEvoked,permTestResults,...
                ~,ROIsInBarrel,~,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_NH{K} );
            traceByStim=traceByStim(2);
            sponTrace=sponTrace(2);
            permTestResults=permTestResults(2);
            
    end
        
    if ~isempty(ROIs_NH)
        sigROIs=ROIs_NH{K};
    else
        sigROIs=find_sigROIs(permTestResults,traceByStim);
    end
    results_field_NH(K) = find_corr_simple( traceByStim,sponTrace,sigROIs,framesEvoked,ROI_positions,mag);
%     [ results_SW_NH(K) ] = find_corrToSWs(sigROIs,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
%     
%     if ~isempty(results_SW_NH(K).ROIdist_inArc)
%        figure(f2)
%         NC_plot=mean(results_SW_NH(K).NC_inArc,2);
%        
%         plot_scatterRLine_color(results_SW_NH(K).ROIdist_inArc,NC_plot,cmap(K,:))
%     end
    
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

% SC_arc_NH=cat(1,results_SW_NH(:).SC_inArc);
% SC_row_NH=cat(1,results_SW_NH(:).SC_inRow);
% NC_arc_NH=cat(1,results_SW_NH(:).NC_inArc);
% NC_arc_NH=mean(NC_arc_NH,2);
% NC_row_NH=cat(1,results_SW_NH(:).NC_inRow);
% NC_row_NH=mean(NC_row_NH,2);
% ROIdist_row_NH=cat(2,results_SW_NH(:).ROIdist_inRow);
% ROIdist_arc_NH=cat(2,results_SW_NH(:).ROIdist_inArc);
% 
% SC_arc_EN=cat(1,results_SW_EN(:).SC_inArc);
% SC_row_EN=cat(1,results_SW_EN(:).SC_inRow);
% NC_arc_EN=cat(1,results_SW_EN(:).NC_inArc);
% NC_arc_EN=mean(NC_arc_EN,2);
% NC_row_EN=cat(1,results_SW_EN(:).NC_inRow);
% NC_row_EN=mean(NC_row_EN,2);
% ROIdist_row_EN=cat(2,results_SW_EN(:).ROIdist_inRow);
% ROIdist_arc_EN=cat(2,results_SW_EN(:).ROIdist_inArc);
% 
% SW_EN.SC_arc=SC_arc_EN;
% SW_EN.SC_row=SC_row_EN;
% SW_EN.NC_arc=NC_arc_EN;
% SW_EN.NC_row=NC_row_EN;
% SW_EN.dist_arc=ROIdist_arc_EN;
% SW_EN.dist_row=ROIdist_row_EN;
% 
% SW_NH.SC_arc=SC_arc_NH;
% SW_NH.SC_row=SC_row_NH;
% SW_NH.NC_arc=NC_arc_NH;
% SW_NH.NC_row=NC_row_NH;
% SW_NH.dist_arc=ROIdist_arc_NH;
% SW_NH.dist_row=ROIdist_row_NH;
SW_NH=[];
SW_EN=[]
%% make plots comparing SC,NC by distance within column pairs vs all

% find NC,SC for ROIs within and across columns
results_NH=[];
for K=1:length(paths_NH)
    switch npSub
        case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_NH{K} );
        case 1
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_NH{K} );
            traceByStim=traceByStim(2);
            sponTrace=sponTrace(2);
            permTestResults=permTestResults(2);
    end
        
    if ~isempty(ROIs_NH)
        sigROIs=ROIs_NH{K};
    else
        sigROIs=find_sigROIs(permTestResults,traceByStim);
    end
    [ results ] = find_corrInBarrelvsAcross( sigROIs,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
    if ~isempty(results)
        results_NH=[results_NH,results];
    end
end

results_EN=[];
for K=1:length(paths_EN)
    switch npSub
        case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_EN{K} );
        case 1
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_EN{K} );
            traceByStim=traceByStim(2);
            sponTrace=sponTrace(2);
            permTestResults=permTestResults(2);
    end
       
    if ~isempty(ROIs_EN)
        sigROIs=ROIs_EN{K};
    else
        sigROIs=find_sigROIs(permTestResults,traceByStim);
    end
    
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