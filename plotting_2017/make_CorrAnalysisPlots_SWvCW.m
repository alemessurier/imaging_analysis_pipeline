function  make_CorrAnalysisPlots_SWvCW( paths_NH,paths_EN,npSub )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% %% find NC,SC for ROIs with shared CW - enriched
% results_SW_EN=[];
% results_CW_EN=[];
% for K=1:length(paths_EN)
%     switch npSub
%         case 0
%             [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
%                 dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_EN{K} );
%         case 1
%             [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
%                 dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_EN{K} );
%     end
%     [ results_SW ] = find_corrInBarrel_byBW(permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag,'SW');
%     if ~isempty(results_SW(1).NC)
%         results_SW_EN=[results_SW_EN,results_SW];
%     end
%     [ results_CW ] = find_corrInBarrel_byBW(permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag,'CW');
%     if ~isempty(results_CW(1).NC)
%         results_CW_EN=[results_CW_EN,results_CW];
%     end
% end
% 
% %% find NC,SC for ROIs with shared CW - NH
% results_CW_NH=[];
% results_SW_NH=[];
% for K=1:length(paths_NH)
%     switch npSub
%          case 0
%             [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
%                 dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_NH{K} );
%         case 1
%             [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
%                 dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_NH{K} );
%     end
%    [ results_SW ] = find_corrInBarrel_byBW(permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag,'SW');
%     if ~isempty(results_SW(1).NC)
%         results_SW_NH=[results_SW_NH,results_SW];
%     end
%     [ results_CW ] = find_corrInBarrel_byBW(permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag,'CW');
%     if ~isempty(results_CW(1).NC)
%         results_CW_NH=[results_CW_NH,results_CW];
%     end
% end
% 
% 
% 
%  %% make plots of NC and SC by distance, within column, CW-tuned
% 
% NC_NH_CW_inCol=cat(1,results_CW_NH(:).NC);
% NC_NH_CW_inCol=mean(NC_NH_CW_inCol,2);
% shuffNC_NH_CW_inCol=cat(1,results_CW_NH(:).shuffNC);
% shuffNC_NH_CW_inCol=mean(shuffNC_NH_CW_inCol,2);
% 
% ROIdist_NH_CW_inCol=cat(2,results_CW_NH(:).ROIdistance);
% SC_NH_CW_inCol=cat(1,results_CW_NH(:).SC);
% shuffSC_NH_CW_inCol=cat(2,results_CW_NH(:).shuffSC)';
% 
% [numBins_NH_CW_inCol,dist_NH_CW_inCol,meanNC_NH_CW_inCol,semNC_NH_CW_inCol,meanSC_NH_CW_inCol,...
%     semSC_NH_CW_inCol]=binCorrByDist(NC_NH_CW_inCol,shuffNC_NH_CW_inCol,SC_NH_CW_inCol,shuffSC_NH_CW_inCol,ROIdist_NH_CW_inCol);
% 
%    
%         NC_EN_CW_inCol=cat(1,results_CW_EN(:).NC);
% NC_EN_CW_inCol=mean(NC_EN_CW_inCol,2);
% shuffNC_EN_CW_inCol=cat(1,results_CW_EN(:).shuffNC);
% shuffNC_EN_CW_inCol=mean(shuffNC_EN_CW_inCol,2);
% 
% ROIdist_EN_CW_inCol=cat(2,results_CW_EN(:).ROIdistance);
% SC_EN_CW_inCol=cat(1,results_CW_EN(:).SC);
% shuffSC_EN_CW_inCol=cat(2,results_CW_EN(:).shuffSC)';
% 
% 
% [numBins_EN_CW_inCol,dist_EN_CW_inCol,meanNC_EN_CW_inCol,semNC_EN_CW_inCol,meanSC_EN_CW_inCol,...
%     semSC_EN_CW_inCol]=binCorrByDist(NC_EN_CW_inCol,shuffNC_EN_CW_inCol,SC_EN_CW_inCol,shuffSC_EN_CW_inCol,ROIdist_EN_CW_inCol);
% 
% 
% 
% % make CDF plots
% compare_cumDists(NC_NH_CW_inCol,NC_EN_CW_inCol,'NC in column, CW-tuned')
% compare_cumDists(SC_NH_CW_inCol,SC_EN_CW_inCol,'SC in column, CW-tuned')
% 
% 
% 
% 
% %  make plots for within column data
% 
% figure;
% ax=make_distVCorr_plot(dist_NH_CW_inCol,numBins_NH_CW_inCol,dist_EN_CW_inCol,numBins_EN_CW_inCol,...
%             meanNC_NH_CW_inCol,meanNC_EN_CW_inCol,semNC_NH_CW_inCol,semNC_EN_CW_inCol);
%         
% ylabel('mean NC')
% title('within column pairs, CW-tuned')
% hold on
% 
% 
% figure;
% make_distVCorr_plot(dist_NH_CW_inCol,numBins_NH_CW_inCol,dist_EN_CW_inCol,numBins_EN_CW_inCol,...
%             meanSC_NH_CW_inCol,meanSC_EN_CW_inCol,semSC_NH_CW_inCol,semSC_EN_CW_inCol)
% ylabel('mean SC')
% title('within column pairs, CW-tuned')
% 
%  %% make plots of NC and SC by distance, within column, SW-tuned
% 
% NC_NH_SW_inCol=cat(1,results_SW_NH(:).NC);
% NC_NH_SW_inCol=mean(NC_NH_SW_inCol,2);
% shuffNC_NH_SW_inCol=cat(1,results_SW_NH(:).shuffNC);
% shuffNC_NH_SW_inCol=mean(shuffNC_NH_SW_inCol,2);
% 
% ROIdist_NH_SW_inCol=cat(2,results_SW_NH(:).ROIdistance);
% SC_NH_SW_inCol=cat(1,results_SW_NH(:).SC);
% shuffSC_NH_SW_inCol=cat(2,results_SW_NH(:).shuffSC)';
% 
% [numBins_NH_SW_inCol,dist_NH_SW_inCol,meanNC_NH_SW_inCol,semNC_NH_SW_inCol,meanSC_NH_SW_inCol,...
%     semSC_NH_SW_inCol]=binCorrByDist(NC_NH_SW_inCol,shuffNC_NH_SW_inCol,SC_NH_SW_inCol,shuffSC_NH_SW_inCol,ROIdist_NH_SW_inCol);
% 
%    
%         NC_EN_SW_inCol=cat(1,results_SW_EN(:).NC);
% NC_EN_SW_inCol=mean(NC_EN_SW_inCol,2);
% shuffNC_EN_SW_inCol=cat(1,results_SW_EN(:).shuffNC);
% shuffNC_EN_SW_inCol=mean(shuffNC_EN_SW_inCol,2);
% 
% ROIdist_EN_SW_inCol=cat(2,results_SW_EN(:).ROIdistance);
% SC_EN_SW_inCol=cat(1,results_SW_EN(:).SC);
% shuffSC_EN_SW_inCol=cat(2,results_SW_EN(:).shuffSC)';
% 
% 
% [numBins_EN_SW_inCol,dist_EN_SW_inCol,meanNC_EN_SW_inCol,semNC_EN_SW_inCol,meanSC_EN_SW_inCol,...
%     semSC_EN_SW_inCol]=binCorrByDist(NC_EN_SW_inCol,shuffNC_EN_SW_inCol,SC_EN_SW_inCol,shuffSC_EN_SW_inCol,ROIdist_EN_SW_inCol);
% 
% 
% 
% % make CDF plots
% compare_cumDists(NC_NH_SW_inCol,NC_EN_SW_inCol,'NC in column, SW-tuned')
% compare_cumDists(SC_NH_SW_inCol,SC_EN_SW_inCol,'SC in column, SW-tuned')
% 
% figure; hold on
% plot_4cdfs(NC_NH_SW_inCol,NC_EN_SW_inCol,NC_NH_CW_inCol,NC_EN_CW_inCol)
% title('within column pairs')
% xlabel('meanNC')
% ylabel('')
% legend('NH,SW-tuned','EN,SW-tuned','NH,CW-tuned','EN,CW-tuned')
% 
% figure; hold on
% plot_4cdfs(SC_NH_SW_inCol,SC_EN_SW_inCol,SC_NH_CW_inCol,SC_EN_CW_inCol)
% title('within column pairs')
% xlabel('mean SC')
% ylabel('')
% legend('NH,SW-tuned','EN,SW-tuned','NH,CW-tuned','EN,CW-tuned')
% 
% 
% %  make plots for within column data
% 
% figure;
% ax=make_distVCorr_plot(dist_NH_SW_inCol,numBins_NH_SW_inCol,dist_EN_SW_inCol,numBins_EN_SW_inCol,...
%             meanNC_NH_SW_inCol,meanNC_EN_SW_inCol,semNC_NH_SW_inCol,semNC_EN_SW_inCol);
%         
% ylabel('mean NC')
% title('within column pairs, SW-tuned')
% hold on
% 
% 
% figure;
% make_distVCorr_plot(dist_NH_SW_inCol,numBins_NH_SW_inCol,dist_EN_SW_inCol,numBins_EN_SW_inCol,...
%             meanSC_NH_SW_inCol,meanSC_EN_SW_inCol,semSC_NH_SW_inCol,semSC_EN_SW_inCol)
% ylabel('mean SC')
% title('within column pairs, SW-tuned')
%% make plots comparing SC,NC by distance within column pairs vs all, CWtuned

% find NC,SC for ROIs within and across columns
results_CW_NH=[];
results_SW_NH=[];
for K=1:length(paths_NH)
    switch npSub
         case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_NH{K} );
        case 1
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_NH{K} );
    end
    [CWtuned,SWtuned] = find_SWtunedCells( ROIsInBarrel,permTestResults,traceByStim,framesEvoked);
     [ results_CW ] = find_corrInBarrelvsAcross( CWtuned,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
     if ~isempty(results_CW)
        results_CW_NH=[results_CW_NH,results_CW];
     end
     
     [ results_SW ] = find_corrInBarrelvsAcross( SWtuned,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
     if ~isempty(results_SW)
        results_SW_NH=[results_SW_NH,results_SW];
     end
end

results_CW_EN=[];
results_SW_EN=[];

for K=1:length(paths_EN)
    switch npSub
         case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_EN{K} );
        case 1
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_EN{K} );
    end
    [CWtuned,SWtuned] = find_SWtunedCells( ROIsInBarrel,permTestResults,traceByStim,framesEvoked);
     [ results_CW ] = find_corrInBarrelvsAcross( CWtuned,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
     if ~isempty(results_CW)
        results_CW_EN=[results_CW_EN,results_CW];
     end
     
     [ results_SW ] = find_corrInBarrelvsAcross( SWtuned,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
     if ~isempty(results_SW)
        results_SW_EN=[results_SW_EN,results_SW];
     end
end


%% make plots for CW-tuned ROIs
  NC_EN_CW_inBarr=cat(1,results_CW_EN(:).NC_inBarr);
NC_EN_CW_inBarr=mean(NC_EN_CW_inBarr,2);
NC_EN_CW_xBarr=cat(1,results_CW_EN(:).NC_xBarr);
NC_EN_CW_xBarr=mean(NC_EN_CW_xBarr,2);
SC_EN_CW_inBarr=cat(1,results_CW_EN(:).SC_inBarr);
SC_EN_CW_xBarr=cat(1,results_CW_EN(:).SC_xBarr);
ROIdist_EN_CW_inBarr=cat(2,results_CW_EN(:).ROIdist_inBarr);
ROIdist_EN_CW_xBarr=cat(2,results_CW_EN(:).ROIdist_xBarr);

 NC_NH_CW_inBarr=cat(1,results_CW_NH(:).NC_inBarr);
NC_NH_CW_inBarr=mean(NC_NH_CW_inBarr,2);
NC_NH_CW_xBarr=cat(1,results_CW_NH(:).NC_xBarr);
NC_NH_CW_xBarr=mean(NC_NH_CW_xBarr,2);
SC_NH_CW_inBarr=cat(1,results_CW_NH(:).SC_inBarr);
SC_NH_CW_xBarr=cat(1,results_CW_NH(:).SC_xBarr);
ROIdist_NH_CW_inBarr=cat(2,results_CW_NH(:).ROIdist_inBarr);
ROIdist_NH_CW_xBarr=cat(2,results_CW_NH(:).ROIdist_xBarr);


[numBins_EN_CW_inBarr,meanDist_EN_CW_inBarr,meanNC_EN_CW_inBarr,semNC_EN_CW_inBarr,...
    meanSC_EN_CW_inBarr,semSC_EN_CW_inBarr]=binCorrByDist(NC_EN_CW_inBarr,NC_EN_CW_inBarr,SC_EN_CW_inBarr,SC_EN_CW_inBarr,ROIdist_EN_CW_inBarr);

[numBins_EN_CW_xBarr,meanDist_EN_CW_xBarr,meanNC_EN_CW_xBarr,semNC_EN_CW_xBarr,...
    meanSC_EN_CW_xBarr,semSC_EN_CW_xBarr]=binCorrByDist(NC_EN_CW_xBarr,NC_EN_CW_xBarr,SC_EN_CW_xBarr,SC_EN_CW_xBarr,ROIdist_EN_CW_xBarr);

[numBins_NH_CW_inBarr,meanDist_NH_CW_inBarr,meanNC_NH_CW_inBarr,semNC_NH_CW_inBarr,...
    meanSC_NH_CW_inBarr,semSC_NH_CW_inBarr]=binCorrByDist(NC_NH_CW_inBarr,NC_NH_CW_inBarr,SC_NH_CW_inBarr,SC_NH_CW_inBarr,ROIdist_NH_CW_inBarr);

[numBins_NH_CW_xBarr,meanDist_NH_CW_xBarr,meanNC_NH_CW_xBarr,semNC_NH_CW_xBarr,...
    meanSC_NH_CW_xBarr,semSC_NH_CW_xBarr]=binCorrByDist(NC_NH_CW_xBarr,NC_NH_CW_xBarr,SC_NH_CW_xBarr,SC_NH_CW_xBarr,ROIdist_NH_CW_xBarr);

%%
figure;
NH_NCplot=make_distVCorr_plot(meanDist_NH_CW_inBarr,numBins_NH_CW_inBarr,meanDist_NH_CW_xBarr,numBins_NH_CW_xBarr,...
            meanNC_NH_CW_inBarr,meanNC_NH_CW_xBarr,semNC_NH_CW_inBarr,semNC_NH_CW_xBarr)
        ylabel('mean NC')
        legend('within column pairs','across column pairs')
        title('Control, CW-tuned')
       
figure;        
NH_SCplot=make_distVCorr_plot(meanDist_NH_CW_inBarr,numBins_NH_CW_inBarr,meanDist_NH_CW_xBarr,numBins_NH_CW_xBarr,...
            meanSC_NH_CW_inBarr,meanSC_NH_CW_xBarr,semSC_NH_CW_inBarr,semSC_NH_CW_xBarr)
        legend('within column pairs','across column pairs')
        ylabel('mean SC')
        title('Control,CW-tuned')
        
figure;
EN_NCplot=make_distVCorr_plot(meanDist_EN_CW_inBarr,numBins_EN_CW_inBarr,meanDist_EN_CW_xBarr,numBins_EN_CW_xBarr,...
            meanNC_EN_CW_inBarr,meanNC_EN_CW_xBarr,semNC_EN_CW_inBarr,semNC_EN_CW_xBarr)
        ylabel('mean NC')
        legend('within column pairs','across column pairs')
        title('Enriched, CW-tuned')
        
        figure;
EN_SCplot=make_distVCorr_plot(meanDist_EN_CW_inBarr,numBins_EN_CW_inBarr,meanDist_EN_CW_xBarr,numBins_EN_CW_xBarr,...
            meanSC_EN_CW_inBarr,meanSC_EN_CW_xBarr,semSC_EN_CW_inBarr,semSC_EN_CW_xBarr)
        legend('within column pairs','across column pairs')
        ylabel('mean SC')
        title('Enriched, CW-tuned')

        
%%
figure;
x_NCplot=make_distVCorr_plot(meanDist_NH_CW_xBarr,numBins_NH_CW_xBarr,meanDist_EN_CW_xBarr,numBins_EN_CW_xBarr,...
            meanNC_NH_CW_xBarr,meanNC_EN_CW_xBarr,semNC_NH_CW_xBarr,semNC_EN_CW_xBarr)
        ylabel('mean NC')
        legend('NH','EN')
        title('across column pairs, CW-tuned')
       
figure;
x_SCplot=make_distVCorr_plot(meanDist_NH_CW_xBarr,numBins_NH_CW_xBarr,meanDist_EN_CW_xBarr,numBins_EN_CW_xBarr,...
            meanSC_NH_CW_xBarr,meanSC_EN_CW_xBarr,semSC_NH_CW_xBarr,semSC_EN_CW_xBarr)
        ylabel('mean SC')
        legend('NH','EN')
        title('across column pairs, CW-tuned')
        
figure;
in_NCplot=make_distVCorr_plot(meanDist_NH_CW_inBarr,numBins_NH_CW_inBarr,meanDist_EN_CW_inBarr,numBins_EN_CW_inBarr,...
            meanNC_NH_CW_inBarr,meanNC_EN_CW_inBarr,semNC_NH_CW_inBarr,semNC_EN_CW_inBarr)
        ylabel('mean NC')
        legend('NH','EN')
        title('within column pairs, CW-tuned')
       
figure;
in_SCplot=make_distVCorr_plot(meanDist_NH_CW_inBarr,numBins_NH_CW_inBarr,meanDist_EN_CW_inBarr,numBins_EN_CW_inBarr,...
            meanSC_NH_CW_inBarr,meanSC_EN_CW_inBarr,semSC_NH_CW_inBarr,semSC_EN_CW_inBarr)
        ylabel('mean SC')
        legend('NH','EN')
        title('within column pairs, CW-tuned')
        



%% make plots for SW-tuned ROIs
  NC_EN_SW_inBarr=cat(1,results_SW_EN(:).NC_inBarr);
NC_EN_SW_inBarr=mean(NC_EN_SW_inBarr,2);
NC_EN_SW_xBarr=cat(1,results_SW_EN(:).NC_xBarr);
NC_EN_SW_xBarr=mean(NC_EN_SW_xBarr,2);
SC_EN_SW_inBarr=cat(1,results_SW_EN(:).SC_inBarr);
SC_EN_SW_xBarr=cat(1,results_SW_EN(:).SC_xBarr);
ROIdist_EN_SW_inBarr=cat(2,results_SW_EN(:).ROIdist_inBarr);
ROIdist_EN_SW_xBarr=cat(2,results_SW_EN(:).ROIdist_xBarr);

 NC_NH_SW_inBarr=cat(1,results_SW_NH(:).NC_inBarr);
NC_NH_SW_inBarr=mean(NC_NH_SW_inBarr,2);
NC_NH_SW_xBarr=cat(1,results_SW_NH(:).NC_xBarr);
NC_NH_SW_xBarr=mean(NC_NH_SW_xBarr,2);
SC_NH_SW_inBarr=cat(1,results_SW_NH(:).SC_inBarr);
SC_NH_SW_xBarr=cat(1,results_SW_NH(:).SC_xBarr);
ROIdist_NH_SW_inBarr=cat(2,results_SW_NH(:).ROIdist_inBarr);
ROIdist_NH_SW_xBarr=cat(2,results_SW_NH(:).ROIdist_xBarr);


[numBins_EN_SW_inBarr,meanDist_EN_SW_inBarr,meanNC_EN_SW_inBarr,semNC_EN_SW_inBarr,...
    meanSC_EN_SW_inBarr,semSC_EN_SW_inBarr]=binCorrByDist(NC_EN_SW_inBarr,NC_EN_SW_inBarr,SC_EN_SW_inBarr,SC_EN_SW_inBarr,ROIdist_EN_SW_inBarr);

[numBins_EN_SW_xBarr,meanDist_EN_SW_xBarr,meanNC_EN_SW_xBarr,semNC_EN_SW_xBarr,...
    meanSC_EN_SW_xBarr,semSC_EN_SW_xBarr]=binCorrByDist(NC_EN_SW_xBarr,NC_EN_SW_xBarr,SC_EN_SW_xBarr,SC_EN_SW_xBarr,ROIdist_EN_SW_xBarr);

[numBins_NH_SW_inBarr,meanDist_NH_SW_inBarr,meanNC_NH_SW_inBarr,semNC_NH_SW_inBarr,...
    meanSC_NH_SW_inBarr,semSC_NH_SW_inBarr]=binCorrByDist(NC_NH_SW_inBarr,NC_NH_SW_inBarr,SC_NH_SW_inBarr,SC_NH_SW_inBarr,ROIdist_NH_SW_inBarr);

[numBins_NH_SW_xBarr,meanDist_NH_SW_xBarr,meanNC_NH_SW_xBarr,semNC_NH_SW_xBarr,...
    meanSC_NH_SW_xBarr,semSC_NH_SW_xBarr]=binCorrByDist(NC_NH_SW_xBarr,NC_NH_SW_xBarr,SC_NH_SW_xBarr,SC_NH_SW_xBarr,ROIdist_NH_SW_xBarr);

%%
% figure;
% NH_NCplot=make_distVCorr_plot(meanDist_NH_SW_inBarr,numBins_NH_SW_inBarr,meanDist_NH_SW_xBarr,numBins_NH_SW_xBarr,...
%             meanNC_NH_SW_inBarr,meanNC_NH_SW_xBarr,semNC_NH_SW_inBarr,semNC_NH_SW_xBarr)
%         ylabel('mean NC')
%         legend('within column pairs','across column pairs')
%         title('Control, SW-tuned')
%        
% figure;        
% NH_SCplot=make_distVCorr_plot(meanDist_NH_SW_inBarr,numBins_NH_SW_inBarr,meanDist_NH_SW_xBarr,numBins_NH_SW_xBarr,...
%             meanSC_NH_SW_inBarr,meanSC_NH_SW_xBarr,semSC_NH_SW_inBarr,semSC_NH_SW_xBarr)
%         legend('within column pairs','across column pairs')
%         ylabel('mean SC')
%         title('Control,SW-tuned')
%         
% figure;
% EN_NCplot=make_distVCorr_plot(meanDist_EN_SW_inBarr,numBins_EN_SW_inBarr,meanDist_EN_SW_xBarr,numBins_EN_SW_xBarr,...
%             meanNC_EN_SW_inBarr,meanNC_EN_SW_xBarr,semNC_EN_SW_inBarr,semNC_EN_SW_xBarr)
%         ylabel('mean NC')
%         legend('within column pairs','across column pairs')
%         title('Enriched, SW-tuned')
%         
%         figure;
% EN_SCplot=make_distVCorr_plot(meanDist_EN_SW_inBarr,numBins_EN_SW_inBarr,meanDist_EN_SW_xBarr,numBins_EN_SW_xBarr,...
%             meanSC_EN_SW_inBarr,meanSC_EN_SW_xBarr,semSC_EN_SW_inBarr,semSC_EN_SW_xBarr)
%         legend('within column pairs','across column pairs')
%         ylabel('mean SC')
%         title('Enriched, SW-tuned')

        
% %%
% figure;
% x_NCplot=make_distVCorr_plot(meanDist_NH_SW_xBarr,numBins_NH_SW_xBarr,meanDist_EN_SW_xBarr,numBins_EN_SW_xBarr,...
%             meanNC_NH_SW_xBarr,meanNC_EN_SW_xBarr,semNC_NH_SW_xBarr,semNC_EN_SW_xBarr)
%         ylabel('mean NC')
%         legend('NH','EN')
%         title('across column pairs, SW-tuned')
%        
% figure;
% x_SCplot=make_distVCorr_plot(meanDist_NH_SW_xBarr,numBins_NH_SW_xBarr,meanDist_EN_SW_xBarr,numBins_EN_SW_xBarr,...
%             meanSC_NH_SW_xBarr,meanSC_EN_SW_xBarr,semSC_NH_SW_xBarr,semSC_EN_SW_xBarr)
%         ylabel('mean SC')
%         legend('NH','EN')
%         title('across column pairs, SW-tuned')
%         
% figure;
% in_NCplot=make_distVCorr_plot(meanDist_NH_SW_inBarr,numBins_NH_SW_inBarr,meanDist_EN_SW_inBarr,numBins_EN_SW_inBarr,...
%             meanNC_NH_SW_inBarr,meanNC_EN_SW_inBarr,semNC_NH_SW_inBarr,semNC_EN_SW_inBarr)
%         ylabel('mean NC')
%         legend('NH','EN')
%         title('within column pairs, SW-tuned')
%        
% figure;
% in_SCplot=make_distVCorr_plot(meanDist_NH_SW_inBarr,numBins_NH_SW_inBarr,meanDist_EN_SW_inBarr,numBins_EN_SW_inBarr,...
%             meanSC_NH_SW_inBarr,meanSC_EN_SW_inBarr,semSC_NH_SW_inBarr,semSC_EN_SW_inBarr)
%         ylabel('mean SC')
%         legend('NH','EN')
%         title('within column pairs, SW-tuned')
% 
% figure; hold on
% plot_4cdfs(NC_NH_SW_inBarr,NC_EN_SW_inBarr,NC_NH_CW_inBarr,NC_EN_CW_inBarr)
% title('within column pairs')
% xlabel('meanNC')
% ylabel('')
% legend('NH,SW-tuned','EN,SW-tuned','NH,CW-tuned','EN,CW-tuned')
% 
% figure; hold on
% plot_4cdfs(SC_NH_SW_inBarr,SC_EN_SW_inBarr,SC_NH_CW_inBarr,SC_EN_CW_inBarr)
% title('within column pairs')
% xlabel('mean SC')
% ylabel('')
% legend('NH,SW-tuned','EN,SW-tuned','NH,CW-tuned','EN,CW-tuned')
% 
% figure; hold on
% plot_4cdfs(NC_NH_SW_xBarr,NC_EN_SW_xBarr,NC_NH_CW_xBarr,NC_EN_CW_xBarr)
% title('across column pairs')
% xlabel('meanNC')
% ylabel('')
% legend('NH,SW-tuned','EN,SW-tuned','NH,CW-tuned','EN,CW-tuned')
% 
% figure; hold on
% plot_4cdfs(SC_NH_SW_xBarr,SC_EN_SW_xBarr,SC_NH_CW_xBarr,SC_EN_CW_xBarr)
% title('across column pairs')
% xlabel('mean SC')
% ylabel('')
% legend('NH,SW-tuned','EN,SW-tuned','NH,CW-tuned','EN,CW-tuned')

%% plot xbarr, inbarr on same axis

figure; hold on
make_distVCorr_plot4(meanDist_NH_SW_inBarr,meanDist_EN_SW_inBarr,...
            meanNC_NH_SW_inBarr,meanNC_EN_SW_inBarr,semNC_NH_SW_inBarr,semNC_EN_SW_inBarr,...
            meanDist_NH_SW_xBarr,meanDist_EN_SW_xBarr,...
            meanNC_NH_SW_xBarr,meanNC_EN_SW_xBarr,semNC_NH_SW_xBarr,semNC_EN_SW_xBarr)
ylabel('mean NC')
title('SW-tuned')
        legend('NH in col','EN in col','NH x col','EN x col')


figure; hold on
make_distVCorr_plot4(meanDist_NH_SW_inBarr,meanDist_EN_SW_inBarr,...
            meanSC_NH_SW_inBarr,meanSC_EN_SW_inBarr,semSC_NH_SW_inBarr,semSC_EN_SW_inBarr,...
            meanDist_NH_SW_xBarr,meanDist_EN_SW_xBarr,...
            meanSC_NH_SW_xBarr,meanSC_EN_SW_xBarr,semSC_NH_SW_xBarr,semSC_EN_SW_xBarr)
ylabel('mean SC')
title('SW-tuned')
        legend('NH in col','EN in col','NH x col','EN x col')


figure; hold on
make_distVCorr_plot4(meanDist_NH_CW_inBarr,meanDist_EN_CW_inBarr,...
            meanSC_NH_CW_inBarr,meanSC_EN_CW_inBarr,semSC_NH_CW_inBarr,semSC_EN_CW_inBarr,...
            meanDist_NH_CW_xBarr,meanDist_EN_CW_xBarr,...
            meanSC_NH_CW_xBarr,meanSC_EN_CW_xBarr,semSC_NH_CW_xBarr,semSC_EN_CW_xBarr)
ylabel('mean SC')
title('CW-tuned')
        legend('NH in col','EN in col','NH x col','EN x col')


figure; hold on
make_distVCorr_plot4(meanDist_NH_CW_inBarr,meanDist_EN_CW_inBarr,...
            meanNC_NH_CW_inBarr,meanNC_EN_CW_inBarr,semNC_NH_CW_inBarr,semNC_EN_CW_inBarr,...
            meanDist_NH_CW_xBarr,meanDist_EN_CW_xBarr,...
            meanNC_NH_CW_xBarr,meanNC_EN_CW_xBarr,semNC_NH_CW_xBarr,semNC_EN_CW_xBarr)
ylabel('mean NC')
title('CW-tuned')
        legend('NH in col','EN in col','NH x col','EN x col')
        
%% Plot CW, SW tuned on same axis

figure; hold on
make_distVCorr_plot4(meanDist_NH_CW_inBarr,meanDist_EN_CW_inBarr,...
            meanNC_NH_CW_inBarr,meanNC_EN_CW_inBarr,semNC_NH_CW_inBarr,semNC_EN_CW_inBarr,...
            meanDist_NH_SW_inBarr,meanDist_EN_SW_inBarr,...
            meanNC_NH_SW_inBarr,meanNC_EN_SW_inBarr,semNC_NH_SW_inBarr,semNC_EN_SW_inBarr)
ylabel('mean NC')
title('Within Barrel pairs')
        legend('NH CW','EN CW','NH SW','EN SW')
        
        
figure; hold on
make_distVCorr_plot4(meanDist_NH_CW_inBarr,meanDist_EN_CW_inBarr,...
            meanSC_NH_CW_inBarr,meanSC_EN_CW_inBarr,semSC_NH_CW_inBarr,semSC_EN_CW_inBarr,...
            meanDist_NH_SW_inBarr,meanDist_EN_SW_inBarr,...
            meanSC_NH_SW_inBarr,meanSC_EN_SW_inBarr,semSC_NH_SW_inBarr,semSC_EN_SW_inBarr)
ylabel('mean SC')
title('Within Barrel pairs')
        legend('NH CW','EN CW','NH SW','EN SW')

        
        figure; hold on
make_distVCorr_plot4(meanDist_NH_CW_xBarr,meanDist_EN_CW_xBarr,...
            meanNC_NH_CW_xBarr,meanNC_EN_CW_xBarr,semNC_NH_CW_xBarr,semNC_EN_CW_xBarr,...
            meanDist_NH_SW_xBarr,meanDist_EN_SW_xBarr,...
            meanNC_NH_SW_xBarr,meanNC_EN_SW_xBarr,semNC_NH_SW_xBarr,semNC_EN_SW_xBarr)
ylabel('mean NC')
title('Across Barrel pairs')
        legend('NH CW','EN CW','NH SW','EN SW')
        
        
figure; hold on
make_distVCorr_plot4(meanDist_NH_CW_xBarr,meanDist_EN_CW_xBarr,...
            meanSC_NH_CW_xBarr,meanSC_EN_CW_xBarr,semSC_NH_CW_xBarr,semSC_EN_CW_xBarr,...
            meanDist_NH_SW_xBarr,meanDist_EN_SW_xBarr,...
            meanSC_NH_SW_xBarr,meanSC_EN_SW_xBarr,semSC_NH_SW_xBarr,semSC_EN_SW_xBarr)
ylabel('mean SC')
title('Across Barrel pairs')
        legend('NH CW','EN CW','NH SW','EN SW')


%% subfunction definitions  

    function [numBins,meanDist,meanNC,semNC,meanSC,semSC,shuff_meanNC,...
            shuff_semNC,shuff_meanSC,shuff_semSC]=binCorrByDist(NC,shuffNC,SC,shuffSC,ROIdist)
        
        
        [ROIdist,sortInds]=sort(ROIdist,'ascend');
        NC=NC(sortInds);
        SC=SC(sortInds);
        shuffNC=shuffNC(sortInds);
        shuffSC=shuffSC(sortInds);
        
        shuffSC=shuffSC(sortInds);
        [numBins,binEdges,binInds]=histcounts(ROIdist,20);
        
        for i=1:20
            meanDist(i)=mean(ROIdist(binInds==i));
            
            meanNC(i)=mean(NC(binInds==i));
            stdNC(i)=std(NC(binInds==i));
            semNC(i)=stdNC(i)/sqrt(sum(binInds==i));
            
            
            shuff_meanNC(i)=mean(shuffNC(binInds==i));
            shuff_stdNC(i)=std(shuffNC(binInds==i));
            shuff_semNC(i)=shuff_stdNC(i)/sqrt(sum(binInds==i));
            
            meanSC(i)=mean(SC(binInds==i));
            stdSC(i)=std(SC(binInds==i));
            semSC(i)=stdSC(i)/sqrt(sum(binInds==i));
            
            shuff_stdSC(i)=std(shuffSC(binInds==i));
            shuff_semSC(i)=shuff_stdSC(i)/sqrt(sum(binInds==i));
        end
    end

    function ax=make_distVCorr_plot(dist_NH,numBins_NH,dist_EN,numBins_EN,...
            meanCorr_NH,meanCorr_EN,semCorr_NH,semCorr_EN)
       
        
        indsUse_EN=numBins_EN>10;
        indsUse_NH=numBins_NH>10;
       
        numBins_NH=numBins_NH(indsUse_NH);
        dist_NH=dist_NH(indsUse_NH);
        meanCorr_NH=meanCorr_NH(indsUse_NH);
        semCorr_NH=semCorr_NH(indsUse_NH);
  
        numBins_EN=numBins_EN(indsUse_EN);
        dist_EN=dist_EN(indsUse_EN);
        meanCorr_EN=meanCorr_EN(indsUse_EN);
        semCorr_EN=semCorr_EN(indsUse_EN);
  
        ax(2)=axes('XAxisLocation','bottom',...
            'YAxisLocation','right',...
            'Color','none',...
            'Position',[0.2    0.1106    0.70    0.8144]);
        hold on
        plot(dist_NH,numBins_NH,'k:');
        plot(dist_EN,numBins_EN,'r:');
        
        ylabel('# pairs')
        
        
        position=ax(2).Position;
        ax(1) = axes('Position',position,...
            'Color','none');
        ax(1).XLim=ax(2).XLim;
        hold on
        ax(2).XTickLabel=[];
        [hl,hp]=boundedline(dist_NH,meanCorr_NH,semCorr_NH,'ko-',...
            dist_EN,meanCorr_EN,semCorr_EN,'ro-',...
            'alpha');
        linkaxes(ax,'x')
        xlabel('distance between ROIs')
        hp_info=get(hp);
        hp_info(1).Annotation.LegendInformation.IconDisplayStyle='off';
        hp_info(2).Annotation.LegendInformation.IconDisplayStyle='off';
        legend('NH','EN')
    end

function make_distVCorr_plot4(meanDist_NH_inBarr,meanDist_EN_inBarr,...
            mean_NH_inBarr,mean_EN_inBarr,sem_NH_inBarr,sem_EN_inBarr,...
            meanDist_NH_xBarr,meanDist_EN_xBarr,...
            mean_NH_xBarr,mean_EN_xBarr,sem_NH_xBarr,sem_EN_xBarr)
       
        
        indsUse_EN_inBarr=~isnan(meanDist_EN_inBarr);
        indsUse_NH_inBarr=~isnan(meanDist_NH_inBarr);
       indsUse_EN_xBarr=~isnan(meanDist_EN_xBarr);
        indsUse_NH_xBarr=~isnan(meanDist_NH_xBarr);
       
        
        
        meanDist_NH_inBarr=meanDist_NH_inBarr(indsUse_NH_inBarr);
        mean_NH_inBarr=mean_NH_inBarr(indsUse_NH_inBarr);
        sem_NH_inBarr=sem_NH_inBarr(indsUse_NH_inBarr);
  
         meanDist_EN_inBarr=meanDist_EN_inBarr(indsUse_EN_inBarr);
        mean_EN_inBarr=mean_EN_inBarr(indsUse_EN_inBarr);
        sem_EN_inBarr=sem_EN_inBarr(indsUse_EN_inBarr);
        
         meanDist_NH_xBarr=meanDist_NH_xBarr(indsUse_NH_xBarr);
        mean_NH_xBarr=mean_NH_xBarr(indsUse_NH_xBarr);
        sem_NH_xBarr=sem_NH_xBarr(indsUse_NH_xBarr);
  
         meanDist_EN_xBarr=meanDist_EN_xBarr(indsUse_EN_xBarr);
        mean_EN_xBarr=mean_EN_xBarr(indsUse_EN_xBarr);
        sem_EN_xBarr=sem_EN_xBarr(indsUse_EN_xBarr);
     
  
%         ax(2)=axes('XAxisLocation','bottom',...
%             'YAxisLocation','right',...
%             'Color','none',...
%             'Position',[0.2    0.1106    0.70    0.8144]);
%         hold on
%         plot(dist_NH,numBins_NH,'k:');
%         plot(dist_EN,numBins_EN,'r:');
%         
%         ylabel('# pairs')
%         
%         
%         position=ax(2).Position;
%         ax(1) = axes('Position',[0.2    0.1106    0.70    0.8144],...
%             'Color','none');
%    
        hold on
   
        [hl,hp]=boundedline(meanDist_NH_inBarr,mean_NH_inBarr,sem_NH_inBarr,'ko-',...
            meanDist_EN_inBarr,mean_EN_inBarr,sem_EN_inBarr,'ro-',...
            meanDist_NH_xBarr,mean_NH_xBarr,sem_NH_xBarr,'ko:',...
            meanDist_EN_xBarr,mean_EN_xBarr,sem_EN_xBarr,'ro:',...
            'alpha');
        
        xlabel('distance between ROIs')
        hp_info=get(hp);
        hp_info(1).Annotation.LegendInformation.IconDisplayStyle='off';
        hp_info(2).Annotation.LegendInformation.IconDisplayStyle='off';
    end

    function compare_cumDists(NH,EN,titleStr)
        figure; hold on
        p1=cdfplot(NH);
        p2=cdfplot(EN);
        p1.Color='k';
        p2.Color='r';
        p1.LineWidth=2.5;
        p2.LineWidth=2.5;
        xlabel('corr coeff')
        ylabel('');
        pval=permutationTest(NH,EN,10000);
        legend('NH','EN')
        title([titleStr, 'p=',num2str(pval)])
    end
end

