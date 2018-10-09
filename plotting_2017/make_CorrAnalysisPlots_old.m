function  [pvals_NC,pvals_SC]=make_CorrAnalysisPlots_old( paths_NH,paths_EN,npSub )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% find NC,SC for ROIs with shared CW - enriched
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
    sigROIs=find_sigROIs(permTestResults,traceByStim);
    [ results ] = find_corrToSWs(sigROIs,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag)
    [ results ] = find_corrInBarrel(permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag );
    if ~isempty(results(1).NC)
        results_EN=[results_EN,results];
    end
end

%% find NC,SC for ROIs with shared CW - NH
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
     [ results ] = find_corrInBarrel(permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag );
      if ~isempty(results(1).NC)
        results_EN=[results_EN,results];
    end
    results_NH=[results_NH,results];
end


%% same BW

results_field_EN=[];

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
    [ results ] = find_corr_sameBW( traceByStim,sponTrace,permTestResults,framesEvoked,ROI_positions,mag);
    results_field_EN=[results_field_EN,results];
end

%% same BW

results_field_NH=[];

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
    [ results ] = find_corr_sameBW( traceByStim,sponTrace,permTestResults,framesEvoked,ROI_positions,mag);
    results_field_NH=[results_field_NH,results];
end

%% plot SC vs NC within barrels, by field

% EN
for K=1:length(results_field_EN)
    NCs=mean(results_field_EN(K).NC,2);
    SCs=results_field_EN(K).SC;
    [r2_EN(K),p_EN(K),slope_EN(K)]=plot_scatterRLine(SCs,NCs);
    xlabel('EN, SC')
    ylabel('EN, NC')
end

%NH
for K=1:length(results_field_NH)
    NCs=mean(results_field_NH(K).NC,2);
    SCs=results_field_NH(K).SC;
    [r2_NH(K),p_NH(K),slope_NH(K)]=plot_scatterRLine(SCs,NCs);
    xlabel('NH, SC')
    ylabel('NH, NC')
end

[ pval] = notBoxPlot_ENvNH( slope_EN,slope_NH );

%%
NC_shared_NH=cat(1,results_field_NH(:).NC_shared);
NC_shared_NH=mean(NC_shared_NH,2);

NC_NH_all=cat(1,results_field_NH(:).NC);
NC_NH_all=mean(NC_NH_all,2);

NC_shared_EN=cat(1,results_field_EN(:).NC_shared);
NC_shared_EN=mean(NC_shared_EN,2);

NC_EN_all=cat(1,results_field_EN(:).NC);
NC_EN_all=mean(NC_EN_all,2);

figure; hold on
plot_4cdfs(NC_shared_NH,NC_shared_EN,NC_NH_all,NC_EN_all)
legend('BWshared NH','BWshared EN','NH all','ENall')
title('NC for shared BW pairs vs. all pairs')
xlabel('mean noise corr across stims')
ylabel('')

%%
SC_shared_NH=arrayfun(@(x)x.SC(x.samePWinds),results_field_NH,'Uni',0);
SC_shared_NH=cat(1,SC_shared_NH{:});
SC_NH_all=cat(1,results_field_NH(:).SC);


SC_shared_EN=arrayfun(@(x)x.SC(x.samePWinds),results_field_EN,'Uni',0);
SC_shared_EN=cat(1,SC_shared_EN{:});
SC_EN_all=cat(1,results_field_EN(:).SC);

figure; hold on
plot_4cdfs(SC_shared_NH,SC_shared_EN,SC_NH_all,SC_EN_all)
legend('BWshared NH','BWshared EN','NH all','ENall')
title('SC for shared BW pairs vs. all pairs')
xlabel('signal corr')
ylabel('')

%% make plots of NC and SC by distance, whole field

NC_NH_all=cat(1,results_field_NH(:).NC);
NC_NH_all=mean(NC_NH_all,2);
shuffNC_NH_all=cat(1,results_field_NH(:).shuffNC);
shuffNC_NH_all=mean(shuffNC_NH_all,2);

ROIdist_NH_all=cat(2,results_field_NH(:).ROIdistance);
SC_NH_all=cat(1,results_field_NH(:).SC);
shuffSC_NH_all=cat(2,results_field_NH(:).shuffSC)';

[numBins_NH,binEdges_NH,dist_NH,meanNC_NH,semNC_NH,meanSC_NH,semSC_NH,allNC_NH,allSC_NH]=...
    binCorrByDist(NC_NH_all,shuffNC_NH_all,SC_NH_all,shuffSC_NH_all,ROIdist_NH_all,20);
        
     
        
%%


NC_EN_all=cat(1,results_field_EN(:).NC);
NC_EN_all=mean(NC_EN_all,2);
shuffNC_EN_all=cat(1,results_field_EN(:).shuffNC);
shuffNC_EN_all=mean(shuffNC_EN_all,2);

ROIdist_EN_all=cat(2,results_field_EN(:).ROIdistance);
SC_EN_all=cat(1,results_field_EN(:).SC);
shuffSC_EN_all=cat(2,results_field_EN(:).shuffSC)';

[numBins_EN,~,dist_EN,meanNC_EN,semNC_EN,meanSC_EN,semSC_EN,allNC_EN,allSC_EN]=...
    binCorrByDist(NC_EN_all,shuffNC_EN_all,SC_EN_all,shuffSC_EN_all,ROIdist_EN_all,binEdges_NH);
 
for i=1:length(allNC_NH)
    tmp_p(i)=permutationTest_mean(allNC_NH{i},allNC_EN{i},10000);
    tmp_p2(i)=permutationTest_mean(allSC_NH{i},allSC_EN{i},10000);
end
pvals_NC.field=tmp_p;
pvals_SC.field=tmp_p2;

% make CDF plots
compare_cumDists(NC_NH_all,NC_EN_all,'NC whole field,')
compare_cumDists(SC_NH_all,SC_EN_all,'SC whole field,')
% compare_cumDists(NC_NH_all,shuffNC_NH_all,'NC whole field,')
% legend('NH','shuffled')
% compare_cumDists(NC_EN_all,shuffNC_EN_all,'NC whole field,')
% legend('EN','shuffled')
% compare_cumDists(SC_NH_all,shuffSC_NH_all,'SC whole field,')
% legend('NH','shuffled')
% compare_cumDists(SC_EN_all,shuffSC_EN_all,'SC whole field,')
% legend('EN','shuffled')


%%  make plots

figure;
make_distVCorr_plot(dist_NH,numBins_NH,dist_EN,numBins_EN,...
            meanNC_NH,meanNC_EN,semNC_NH,semNC_EN)
        
ylabel('mean NC')

figure;
make_distVCorr_plot(dist_NH,numBins_NH,dist_EN,numBins_EN,...
            meanSC_NH,meanSC_EN,semSC_NH,semSC_EN)
ylabel('mean SC')

 %% make plots of NC and SC by distance, within column

NC_NH_CW=cat(1,results_NH(:).NC);
NC_NH_CW=mean(NC_NH_CW,2);
shuffNC_NH_CW=cat(1,results_NH(:).shuffNC);
shuffNC_NH_CW=mean(shuffNC_NH_CW,2);

ROIdist_NH_CW=cat(2,results_NH(:).ROIdistance);
SC_NH_CW=cat(1,results_NH(:).SC);
shuffSC_NH_CW=cat(2,results_NH(:).shuffSC)';
SC_NH_SW=cat(1,results_NH(:).SC_SW);
shuffSC_NH_SW=cat(2,results_NH(:).shuffSC_SW)';
NC_NH_CW_only=cat(1,results_NH(:).NC_CW);
shuffNC_NH_CW_only=cat(1,results_NH(:).shuffNC_CW);

[numBins_NH_CW,binEdges_NH,dist_NH_CW,meanNC_NH_CW,semNC_NH_CW,meanSC_NH_CW,...
    semSC_NH_CW,colNC_NH,colSC_NH]=binCorrByDist(NC_NH_CW,shuffNC_NH_CW,SC_NH_CW,shuffSC_NH_CW,ROIdist_NH_CW,20);

[~,~,~,meanNC_NH_CW_only,semNC_NH_CW_only]=binCorrByDist(NC_NH_CW_only,shuffNC_NH_CW_only,SC_NH_CW,shuffSC_NH_CW,ROIdist_NH_CW,20);

[~,~,~,~,~,meanSC_NH_SW,semSC_NH_SW]=binCorrByDist(NC_NH_CW,shuffNC_NH_CW,SC_NH_SW,shuffSC_NH_SW,ROIdist_NH_CW,20);
        %%
        NC_EN_CW=cat(1,results_EN(:).NC);
NC_EN_CW=mean(NC_EN_CW,2);
shuffNC_EN_CW=cat(1,results_EN(:).shuffNC);
shuffNC_EN_CW=mean(shuffNC_EN_CW,2);

ROIdist_EN_CW=cat(2,results_EN(:).ROIdistance);
SC_EN_CW=cat(1,results_EN(:).SC);
shuffSC_EN_CW=cat(2,results_EN(:).shuffSC)';
SC_EN_SW=cat(1,results_EN(:).SC_SW);
shuffSC_EN_SW=cat(2,results_EN(:).shuffSC_SW)';
NC_EN_CW_only=cat(1,results_EN(:).NC_CW);
shuffNC_EN_CW_only=cat(1,results_EN(:).shuffNC_CW);


[numBins_EN_CW,~,dist_EN_CW,meanNC_EN_CW,semNC_EN_CW,meanSC_EN_CW,...
    semSC_EN_CW,colNC_EN,colSC_EN]=binCorrByDist(NC_EN_CW,shuffNC_EN_CW,SC_EN_CW,shuffSC_EN_CW,ROIdist_EN_CW,binEdges_NH);

[~,~,~,meanNC_EN_CW_only,semNC_EN_CW_only]=binCorrByDist(NC_EN_CW_only,shuffNC_EN_CW_only,SC_EN_CW,shuffSC_EN_CW,ROIdist_EN_CW,20);

[~,~,~,~,~,meanSC_EN_SW,semSC_EN_SW]=binCorrByDist(NC_EN_CW,shuffNC_EN_CW,SC_EN_SW,shuffSC_EN_SW,ROIdist_EN_CW,20);


% make CDF plots
compare_cumDists(NC_NH_CW,NC_EN_CW,'NC in column,')
compare_cumDists(SC_NH_CW,SC_EN_CW,'SC in column,')
% compare_cumDists(SC_NH_SW,SC_EN_SW,'SC in column, surround RF,')

compare_cumDists(NC_NH_CW,NC_NH_all,'NC NH field vs col');
compare_cumDists(NC_EN_CW,NC_EN_all,'NC EN field vs col');

compare_cumDists(SC_NH_CW,SC_NH_all,'SC NH field vs col');
compare_cumDists(SC_EN_CW,SC_EN_all,'SC EN field vs col');

% make CDF plots, within column data and full-field data in same axes
figure; hold on
plot_4cdfs(NC_NH_all,NC_EN_all,NC_NH_CW,NC_EN_CW)
legend('NH,field','EN,field','NH,column','EN,column')
ylabel('')
xlabel('pairwise noise corr')
title('NC')

figure; hold on
plot_4cdfs(SC_NH_all,SC_EN_all,SC_NH_CW,SC_EN_CW)
legend('NH,field','EN,field','NH,column','EN,column')
ylabel('')
xlabel('pairwise signal corr')
title('SC')
% compare_cumDists(NC_NH_CW,shuffNC_NH_CW,'NC in column,')
% legend('NH','shuffled')
% compare_cumDists(NC_EN_CW,shuffNC_EN_CW,'NC in column,')
% legend('EN','shuffled')
% compare_cumDists(SC_NH_CW,shuffSC_NH_CW,'SC in column,')
% legend('NH','shuffled')
% compare_cumDists(SC_EN_CW,shuffSC_EN_CW,'SC in column,')
% legend('EN','shuffled')
% compare_cumDists(SC_NH_SW,shuffSC_NH_SW,'SC in column, surround RF,')
% legend('NH','shuffled')
% compare_cumDists(SC_EN_SW,shuffSC_EN_SW,'SC in column,surround RF')
% legend('EN','shuffled')

%%  make plots for within column data

figure;
ax=make_distVCorr_plot(dist_NH_CW,numBins_NH_CW,dist_EN_CW,numBins_EN_CW,...
            meanNC_NH_CW,meanNC_EN_CW,semNC_NH_CW,semNC_EN_CW);
        
ylabel('mean NC')
title('within column pairs')
hold on


tmp_p=zeros(1,length(colNC_NH));
tmp_p2=zeros(1,length(colNC_NH))
for i=1:length(colNC_NH)
    if numBins_NH_CW(i)>2 && numBins_EN_CW(i)>2
        tmp_p(i)=permutationTest_mean(colNC_NH{i},colNC_EN{i},10000);
        tmp_p2(i)=permutationTest_mean(colSC_NH{i},colSC_EN{i},10000);
    else
        tmp_p(i)=nan;
        tmp_p2(i)=nan;
    end
end
pvals_NC.col=tmp_p;
pvals_SC.col=tmp_p2;
% % make rotated histogram of NC
% 
% tmp2=axes('Position',[0.05    0.1191    0.0653    0.8059],...
%     'Units','normalized');
%  binMin=ax(1).YLim(1);
%  binMax=ax(1).YLim(2);
%  binEdges=binMin:(binMax-binMin)/20:binMax;
%  
%  histNC_NH=histc(NC_NH_CW,binEdges);
%  histNC_NH=histNC_NH/sum(histNC_NH);
%  histNC_EN=histc(NC_EN_CW,binEdges);
% histNC_EN=histNC_EN/sum(histNC_EN);
% plot(histNC_NH,binEdges,'k-','LineWidth',1.5);
% hold on
% plot(histNC_EN,binEdges,'r-','LineWidth',1.5);
% %  handlebar_NH=barh(binEdges,histNC_NH);
% % handlebar_NH.FaceColor='none';
% % handlebar_NH.EdgeColor='k';
% % handlebar_NH.LineWidth=1.5;
% % hold on
% % handlebar_EN=barh(binEdges,histNC_EN);
% % handlebar_EN.FaceColor='none';
% % handlebar_EN.EdgeColor='r';
% % handlebar_EN.LineWidth=1.5;
% tmp2=gca;
% tmp2.XDir='reverse';
% tmp2.YLim=ax(1).YLim;
% tmp2.YTick=[];
% % tmp2.Position=[0.05    0.1191    0.0653    0.8059];
% tmp2.XLim=[0 max([histNC_NH;histNC_EN])];
% tmp2.XTick=[0 max([histNC_NH;histNC_EN])];
% tmp2.XTickLabel=[0 max([histNC_NH;histNC_EN])];
% tmp2.YAxisLocation='right';
% tmp2.Box='off';
% tmp.Position=[0.2    0.1106    0.70    0.8144];


figure;
make_distVCorr_plot(dist_NH_CW,numBins_NH_CW,dist_EN_CW,numBins_EN_CW,...
            meanSC_NH_CW,meanSC_EN_CW,semSC_NH_CW,semSC_EN_CW)
ylabel('mean SC')
title('within column pairs')

figure;
make_distVCorr_plot(dist_NH_CW,numBins_NH_CW,dist_EN_CW,numBins_EN_CW,...
            meanNC_NH_CW_only,meanNC_EN_CW_only,semNC_NH_CW_only,semNC_EN_CW_only)
        
ylabel('NC during CW stim')
title('within column pairs')

figure;
make_distVCorr_plot(dist_NH_CW,numBins_NH_CW,dist_EN_CW,numBins_EN_CW,...
            meanSC_NH_SW,meanSC_EN_SW,semSC_NH_SW,semSC_EN_SW)
ylabel('mean SC')
title('in column pairs, surround RF')

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
     [ results ] = find_corrInBarrelvsAcross( permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
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
    [ results ] = find_corrInBarrelvsAcross( permTestResults,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
    if ~isempty(results)
        results_EN=[results_EN,results];
    end
end





%%
  NC_EN_inBarr=cat(1,results_EN(:).NC_inBarr);
NC_EN_inBarr=mean(NC_EN_inBarr,2);
NC_EN_xBarr=cat(1,results_EN(:).NC_xBarr);
NC_EN_xBarr=mean(NC_EN_xBarr,2);
SC_EN_inBarr=cat(1,results_EN(:).SC_inBarr);
SC_EN_xBarr=cat(1,results_EN(:).SC_xBarr);
ROIdist_EN_inBarr=cat(2,results_EN(:).ROIdist_inBarr);
ROIdist_EN_xBarr=cat(2,results_EN(:).ROIdist_xBarr);

 NC_NH_inBarr=cat(1,results_NH(:).NC_inBarr);
NC_NH_inBarr=mean(NC_NH_inBarr,2);
NC_NH_xBarr=cat(1,results_NH(:).NC_xBarr);
NC_NH_xBarr=mean(NC_NH_xBarr,2);
SC_NH_inBarr=cat(1,results_NH(:).SC_inBarr);
SC_NH_xBarr=cat(1,results_NH(:).SC_xBarr);
ROIdist_NH_inBarr=cat(2,results_NH(:).ROIdist_inBarr);
ROIdist_NH_xBarr=cat(2,results_NH(:).ROIdist_xBarr);

[numBins_NH_inBarr,binEdgesInB_NH,meanDist_NH_inBarr,meanNC_NH_inBarr,semNC_NH_inBarr,...
    meanSC_NH_inBarr,semSC_NH_inBarr,inBarrNC_NH,inBarrSC_NH]=binCorrByDist(NC_NH_inBarr,NC_NH_inBarr,SC_NH_inBarr,SC_NH_inBarr,ROIdist_NH_inBarr,20);

[numBins_NH_xBarr,binEdgesXB_NH,meanDist_NH_xBarr,meanNC_NH_xBarr,semNC_NH_xBarr,...
    meanSC_NH_xBarr,semSC_NH_xBarr,xBarrNC_NH,xBarrSC_NH]=binCorrByDist(NC_NH_xBarr,NC_NH_xBarr,SC_NH_xBarr,SC_NH_xBarr,ROIdist_NH_xBarr,20);


[numBins_EN_inBarr,~,meanDist_EN_inBarr,meanNC_EN_inBarr,semNC_EN_inBarr,...
    meanSC_EN_inBarr,semSC_EN_inBarr,inBarrNC_EN,inBarrSC_EN]=binCorrByDist(NC_EN_inBarr,NC_EN_inBarr,SC_EN_inBarr,SC_EN_inBarr,ROIdist_EN_inBarr,binEdgesInB_NH);

[numBins_EN_xBarr,~,meanDist_EN_xBarr,meanNC_EN_xBarr,semNC_EN_xBarr,...
    meanSC_EN_xBarr,semSC_EN_xBarr,xBarrNC_EN,xBarrSC_EN]=binCorrByDist(NC_EN_xBarr,NC_EN_xBarr,SC_EN_xBarr,SC_EN_xBarr,ROIdist_EN_xBarr,binEdgesXB_NH);


% %%
% figure;
% NH_NCplot=make_distVCorr_plot(meanDist_NH_inBarr,numBins_NH_inBarr,meanDist_NH_xBarr,numBins_NH_xBarr,...
%             meanNC_NH_inBarr,meanNC_NH_xBarr,semNC_NH_inBarr,semNC_NH_xBarr)
%         ylabel('mean NC')
%         legend('within column pairs','across column pairs')
%         title('Control')
%        
% figure;        
% NH_SCplot=make_distVCorr_plot(meanDist_NH_inBarr,numBins_NH_inBarr,meanDist_NH_xBarr,numBins_NH_xBarr,...
%             meanSC_NH_inBarr,meanSC_NH_xBarr,semSC_NH_inBarr,semSC_NH_xBarr)
%         legend('within column pairs','across column pairs')
%         ylabel('mean SC')
%         title('Control')
%         
% figure;
% EN_NCplot=make_distVCorr_plot(meanDist_EN_inBarr,numBins_EN_inBarr,meanDist_EN_xBarr,numBins_EN_xBarr,...
%             meanNC_EN_inBarr,meanNC_EN_xBarr,semNC_EN_inBarr,semNC_EN_xBarr)
%         ylabel('mean NC')
%         legend('within column pairs','across column pairs')
%         title('Enriched')
%         
%         figure;
% EN_SCplot=make_distVCorr_plot(meanDist_EN_inBarr,numBins_EN_inBarr,meanDist_EN_xBarr,numBins_EN_xBarr,...
%             meanSC_EN_inBarr,meanSC_EN_xBarr,semSC_EN_inBarr,semSC_EN_xBarr)
%         legend('within column pairs','across column pairs')
%         ylabel('mean SC')
%         title('Enriched')

        
%%
figure;
x_NCplot=make_distVCorr_plot(meanDist_NH_xBarr,numBins_NH_xBarr,meanDist_EN_xBarr,numBins_EN_xBarr,...
            meanNC_NH_xBarr,meanNC_EN_xBarr,semNC_NH_xBarr,semNC_EN_xBarr)
        ylabel('mean NC')
        legend('NH','EN')
        title('across column pairs')
       
figure;
x_SCplot=make_distVCorr_plot(meanDist_NH_xBarr,numBins_NH_xBarr,meanDist_EN_xBarr,numBins_EN_xBarr,...
            meanSC_NH_xBarr,meanSC_EN_xBarr,semSC_NH_xBarr,semSC_EN_xBarr)
        ylabel('mean SC')
        legend('NH','EN')
        title('across column pairs')
        
tmp_p=zeros(1,length(xBarrNC_NH));
tmp_p2=zeros(1,length(xBarrNC_NH))
for i=1:length(xBarrNC_NH)
    if numBins_NH_xBarr(i)>2 && numBins_EN_xBarr(i)>2
        tmp_p(i)=permutationTest_mean(xBarrNC_NH{i},xBarrNC_EN{i},10000);
        tmp_p2(i)=permutationTest_mean(xBarrSC_NH{i},xBarrSC_EN{i},10000);
    else
        tmp_p(i)=nan;
        tmp_p2(i)=nan;
    end
end
pvals_NC.xBarr=tmp_p;
pvals_SC.xBarr=tmp_p2;

%% subfunction definitions  

    function [numBins,binEdges,meanDist,meanNC,semNC,meanSC,semSC,allNC,allSC,shuff_meanNC,...
            shuff_semNC,shuff_meanSC,shuff_semSC]=binCorrByDist(NC,shuffNC,SC,shuffSC,ROIdist,binEdges)
        
        
        [ROIdist,sortInds]=sort(ROIdist,'ascend');
        NC=NC(sortInds);
        SC=SC(sortInds);
        shuffNC=shuffNC(sortInds);
        shuffSC=shuffSC(sortInds);
        
        shuffSC=shuffSC(sortInds);
        [numBins,binEdges,binInds]=histcounts(ROIdist,binEdges);
        
        for i=1:length(numBins)
            meanDist(i)=mean(ROIdist(binInds==i));
            
            allNC{i}=NC(binInds==i);
            meanNC(i)=mean(NC(binInds==i));
            stdNC(i)=std(NC(binInds==i));
            semNC(i)=stdNC(i)/sqrt(sum(binInds==i));
            
            
            shuff_meanNC(i)=mean(shuffNC(binInds==i));
            shuff_stdNC(i)=std(shuffNC(binInds==i));
            shuff_semNC(i)=shuff_stdNC(i)/sqrt(sum(binInds==i));
            
            allSC{i}=SC(binInds==i);
            meanSC(i)=mean(SC(binInds==i));
            stdSC(i)=std(SC(binInds==i));
            semSC(i)=stdSC(i)/sqrt(sum(binInds==i));
            
            shuff_stdSC(i)=std(shuffSC(binInds==i));
            shuff_semSC(i)=shuff_stdSC(i)/sqrt(sum(binInds==i));
        end
    end

    function ax=make_distVCorr_plot(dist_NH,numBins_NH,dist_EN,numBins_EN,...
            meanCorr_NH,meanCorr_EN,semCorr_NH,semCorr_EN)
       
        
        indsUse_EN=numBins_EN>0;
        indsUse_NH=numBins_NH>0;
       
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

