function  [fullfield_NH,fullfield_EN,inRow_NH,inRow_EN,inArc_NH,inArc_EN]=make_CorrAnalysisPlots_surround( paths_NH,paths_EN,npSub,ROIs_NH,ROIs_EN )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if isempty(ROIs_NH)
    for K=1:length(paths_NH)
        switch npSub
            case 0
                [~,traceByStim,~,~,permTestResults] = load_nonNPsub_data( paths_NH{K} );
            case 1
                [~,traceByStim,~,~,permTestResults] = load_nonNPsub_data( paths_NH{K} );
                traceByStim=traceByStim(2);
                sponTrace=sponTrace(2);
                permTestResults=permTestResults(2);
                
        end
        ROIs_NH{K}=find_sigROIs(permTestResults,traceByStim);
    end
end

if isempty(ROIs_EN)
    for K=1:length(paths_EN)
        switch npSub
            case 0
                [~,traceByStim,~,~,permTestResults] = load_nonNPsub_data( paths_EN{K} );
            case 1
                [~,traceByStim,~,~,permTestResults] = load_nonNPsub_data( paths_EN{K} );
                traceByStim=traceByStim(2);
                sponTrace=sponTrace(2);
                permTestResults=permTestResults(2);
                
        end
        ROIs_EN{K}=find_sigROIs(permTestResults,traceByStim);
    end
end

% only use fields that have at least 3 ROIs
useEN=cellfun(@(x)length(x)>2,ROIs_EN);
paths_EN=paths_EN(useEN);
ROIs_EN=ROIs_EN(useEN);

useNH=cellfun(@(x)length(x)>2,ROIs_NH);
paths_NH=paths_NH(useNH);
ROIs_NH=ROIs_NH(useNH);


%% make plots comparing SC,NC by distance within row or within arc (but not within column) pairs

% find NC,SC for ROIs within and across columns
results_NH=[];
for K=1:length(paths_NH)
    switch npSub
        case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsinRowel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_NH{K} );
        case 1
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsinRowel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_NH{K} );
            traceByStim=traceByStim(2);
            sponTrace=sponTrace(2);
            permTestResults=permTestResults(2);
    end
    sigROIs=ROIs_NH{K};
    [ results ] = find_corrToSWs(sigROIs,traceByStim,framesEvoked,ROIsinRowel,sponTrace,ROI_positions,mag);
    if ~isempty(results)
        results_NH=[results_NH,results];
    end
end

results_EN=[];
for K=1:length(paths_EN)
    switch npSub
        case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsinRowel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_EN{K} );
        case 1
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsinRowel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data( paths_EN{K} );
            traceByStim=traceByStim(2);
            sponTrace=sponTrace(2);
            permTestResults=permTestResults(2);
    end
    sigROIs=ROIs_EN{K};
    [ results ] = find_corrToSWs(sigROIs,traceByStim,framesEvoked,ROIsinRowel,sponTrace,ROI_positions,mag);
    if ~isempty(results)
        results_EN=[results_EN,results];
    end
end





%%
NC_EN_inRow=cat(1,results_EN(:).NC_inRow);
inRow_EN.NC=mean(NC_EN_inRow,2);
NC_EN_inArc=cat(1,results_EN(:).NC_inArc);
inArc_EN.NC=mean(NC_EN_inArc,2);
inRow_EN.SC=cat(1,results_EN(:).SC_inRow);
inArc_EN.SC=cat(1,results_EN(:).SC_inArc);
inRow_EN.dist=cat(2,results_EN(:).ROIdist_inRow);
inArc_EN.dist=cat(2,results_EN(:).ROIdist_inArc);



NC_NH_inRow=cat(1,results_NH(:).NC_inRow);
inRow_NH.NC=mean(NC_NH_inRow,2);
NC_NH_inArc=cat(1,results_NH(:).NC_inArc);
inArc_NH.NC=mean(NC_NH_inArc,2);
inRow_NH.SC=cat(1,results_NH(:).SC_inRow);
inArc_NH.SC=cat(1,results_NH(:).SC_inArc);
inRow_NH.dist=cat(2,results_NH(:).ROIdist_inRow);
inArc_NH.dist=cat(2,results_NH(:).ROIdist_inArc);


inArcNC_byROI_NH=cat(1,results_NH(:).inArcNC_byROI);
inRowNC_byROI_NH=cat(1,results_NH(:).inRowNC_byROI);
inColNC_byROI_NH=cat(1,results_NH(:).inColNC_byROI);
inArcNC_byROI_EN=cat(1,results_EN(:).inArcNC_byROI);
inRowNC_byROI_EN=cat(1,results_EN(:).inRowNC_byROI);
inColNC_byROI_EN=cat(1,results_EN(:).inColNC_byROI);

inArcSC_byROI_NH=cat(2,results_NH(:).inArcSC_byROI);
inRowSC_byROI_NH=cat(2,results_NH(:).inRowSC_byROI);
inColSC_byROI_NH=cat(2,results_NH(:).inColSC_byROI);
inArcSC_byROI_EN=cat(2,results_EN(:).inArcSC_byROI);
inRowSC_byROI_EN=cat(2,results_EN(:).inRowSC_byROI);
inColSC_byROI_EN=cat(2,results_EN(:).inColSC_byROI);


%% correlation plots: mean sc in arc v in row by ROI

inds_ArcVrow=~isnan(inArcSC_byROI_NH) & ~isnan(inRowSC_byROI_NH)
ArcVrow_Arc_NH=inArcSC_byROI_NH(inds_ArcVrow);
ArcVrow_Row_NH=inRowSC_byROI_NH(inds_ArcVrow);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVrow_Arc_NH,ArcVrow_Row_NH)
xlabel('SC in arc')
ylabel('SC in ROW (NH)')

inds_ArcVrow=~isnan(inArcSC_byROI_EN) & ~isnan(inRowSC_byROI_EN)
ArcVrow_Arc_EN=inArcSC_byROI_EN(inds_ArcVrow);
ArcVrow_Row_EN=inRowSC_byROI_EN(inds_ArcVrow);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVrow_Arc_EN,ArcVrow_Row_EN)
xlabel('SC in arc');
ylabel('SC in ROW (EN)')
%%
inds_ArcVCol=~isnan(inArcSC_byROI_NH) & ~isnan(inColSC_byROI_NH);
ArcVCol_Arc_NH=inArcSC_byROI_NH(inds_ArcVCol);
ArcVCol_Col_NH=inColSC_byROI_NH(inds_ArcVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVCol_Col_NH,ArcVCol_Arc_NH)
xlabel('SC in column')
ylabel('SC in arc (NH)')

inds_ArcVCol=~isnan(inArcSC_byROI_EN) & ~isnan(inColSC_byROI_EN);
ArcVCol_Arc_EN=inArcSC_byROI_EN(inds_ArcVCol);
ArcVCol_Col_EN=inColSC_byROI_EN(inds_ArcVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVCol_Col_EN,ArcVCol_Arc_EN)
xlabel('SC in column')
ylabel('SC in arc (EN)')
%%
inds_RowVCol=~isnan(inRowSC_byROI_NH) & ~isnan(inColSC_byROI_NH);
RowVCol_Row_NH=inRowSC_byROI_NH(inds_RowVCol);
RowVCol_Col_NH=inColSC_byROI_NH(inds_RowVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(RowVCol_Col_NH,RowVCol_Row_NH)
xlabel('SC in column')
ylabel('SC in ROW (NH)')

inds_RowVCol=~isnan(inRowSC_byROI_EN) & ~isnan(inColSC_byROI_EN);
RowVCol_Row_EN=inRowSC_byROI_EN(inds_RowVCol);
RowVCol_Col_EN=inColSC_byROI_EN(inds_RowVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(RowVCol_Col_EN,RowVCol_Row_EN)
xlabel('SC in column')
ylabel('SC in ROW (EN)')

%% correlation plots: mean nc in arc v in row by ROI

inds_ArcVrow=~isnan(inArcNC_byROI_NH) & ~isnan(inRowNC_byROI_NH);
ArcVrow_Arc_NH=inArcNC_byROI_NH(inds_ArcVrow);
ArcVrow_Row_NH=inRowNC_byROI_NH(inds_ArcVrow);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVrow_Arc_NH,ArcVrow_Row_NH)
xlabel('NC in arc')
ylabel('NC in ROW (NH)')


inds_ArcVrow=~isnan(inArcNC_byROI_EN) & ~isnan(inRowNC_byROI_EN);
ArcVrow_Arc_EN=inArcNC_byROI_EN(inds_ArcVrow);
ArcVrow_Row_EN=inRowNC_byROI_EN(inds_ArcVrow);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVrow_Arc_EN,ArcVrow_Row_EN)
xlabel('NC in arc');
ylabel('NC in ROW (EN)')
%%
inds_ArcVCol=~isnan(inArcNC_byROI_NH) & ~isnan(inColNC_byROI_NH);
ArcVCol_Arc_NH=inArcNC_byROI_NH(inds_ArcVCol);
ArcVCol_Col_NH=inColNC_byROI_NH(inds_ArcVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVCol_Col_NH,ArcVCol_Arc_NH)
xlabel('NC in column')
ylabel('NC in arc (NH)')

inds_ArcVCol=~isnan(inArcNC_byROI_EN) & ~isnan(inColNC_byROI_EN);
ArcVCol_Arc_EN=inArcNC_byROI_EN(inds_ArcVCol);
ArcVCol_Col_EN=inColNC_byROI_EN(inds_ArcVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(ArcVCol_Col_EN,ArcVCol_Arc_EN)
xlabel('NC in column')
ylabel('NC in arc (EN)')
%%
inds_RowVCol=~isnan(inRowNC_byROI_NH) & ~isnan(inColNC_byROI_NH);
RowVCol_Row_NH=inRowNC_byROI_NH(inds_RowVCol);
RowVCol_Col_NH=inColNC_byROI_NH(inds_RowVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(RowVCol_Col_NH,RowVCol_Row_NH)
xlabel('NC in column')
ylabel('NC in ROW (NH)')

inds_RowVCol=~isnan(inRowNC_byROI_EN) & ~isnan(inColNC_byROI_EN);
RowVCol_Row_EN=inRowNC_byROI_EN(inds_RowVCol);
RowVCol_Col_EN=inColNC_byROI_EN(inds_RowVCol);

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(RowVCol_Col_EN,RowVCol_Row_EN)
xlabel('NC in column')
ylabel('NC in ROW (EN)')

%% compare NC vs SC by cell
figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(inArcSC_byROI_EN,inArcNC_byROI_EN)
xlabel('SC in arc'); ylabel('NC in arc (EN)');

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(inArcSC_byROI_NH,inArcNC_byROI_NH)
xlabel('SC in arc'); ylabel('NC in arc (NH)');

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(inRowSC_byROI_EN,inRowNC_byROI_EN)
xlabel('SC in row'); ylabel('NC in row (EN)');

figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(inRowSC_byROI_NH,inRowNC_byROI_NH)
xlabel('SC in row'); ylabel('NC in row (NH)');


figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(inColSC_byROI_EN,inColNC_byROI_EN)
xlabel('SC in col'); ylabel('NC in col (EN)');


figure; hold on
plot([0 1],[0 1],'k:','LineWidth',0.5)
plot_scatterRLine(inColSC_byROI_NH,inColNC_byROI_NH)
xlabel('SC in col'); ylabel('NC in col (NH)');
%% make plots
binEdges=10
[ dist_NH,NC_NH,SEM_NH,num_ROIs_NH,binEdges ] = binVarByDist( inRow_NH.dist,inRow_NH.NC,binEdges );
    [ dist_EN,NC_EN,SEM_EN,num_ROIs_EN,binEdges ] = binVarByDist( inRow_EN.dist,inRow_EN.NC,binEdges );

    totesComboPlot( dist_NH,NC_NH,SEM_NH,num_ROIs_NH,NC_NH,dist_EN,NC_EN,SEM_EN,num_ROIs_EN,NC_EN )
xlabel('distance')
ylabel('NC')
title('In Row')
    %% make plots
binEdges=10
[ dist_NH,NC_NH,SEM_NH,num_ROIs_NH,binEdges ] = binVarByDist( inArc_NH.dist,inArc_NH.NC,binEdges );
    [ dist_EN,NC_EN,SEM_EN,num_ROIs_EN,binEdges ] = binVarByDist( inArc_EN.dist,inArc_EN.NC,binEdges );

    totesComboPlot( dist_NH,NC_NH,SEM_NH,num_ROIs_NH,NC_NH,dist_EN,NC_EN,SEM_EN,num_ROIs_EN,NC_EN )
    xlabel('distance')
ylabel('NC')
title('In Arc')

%% make plots
binEdges=10
[ dist_NH,SC_NH,SEM_NH,num_ROIs_NH,binEdges ] = binVarByDist( inRow_NH.dist,inRow_NH.SC,binEdges );
    [ dist_EN,SC_EN,SEM_EN,num_ROIs_EN,binEdges ] = binVarByDist( inRow_EN.dist,inRow_EN.SC,binEdges );

    totesComboPlot( dist_NH,SC_NH,SEM_NH,num_ROIs_NH,SC_NH,dist_EN,SC_EN,SEM_EN,num_ROIs_EN,SC_EN )
xlabel('distance')
ylabel('SC')
title('In Row')
    %% make plots
binEdges=10
[ dist_NH,NC_NH,SEM_NH,num_ROIs_NH,binEdges ] = binVarByDist( inArc_NH.dist,inArc_NH.SC,binEdges );
    [ dist_EN,NC_EN,SEM_EN,num_ROIs_EN,binEdges ] = binVarByDist( inArc_EN.dist,inArc_EN.SC,binEdges );

    totesComboPlot( dist_NH,NC_NH,SEM_NH,num_ROIs_NH,NC_NH,dist_EN,NC_EN,SEM_EN,num_ROIs_EN,NC_EN )
    xlabel('distance')
ylabel('SC')
title('In Arc')

%%

binEdges=10
[ dist_inRow,NC_inRow,SEM_inRow,num_ROIs_inRow,binEdges ] = binVarByDist( inRow_NH.dist,inRow_NH.NC,binEdges );
    [ dist_inArc,NC_inArc,SEM_inArc,num_ROIs_inArc,binEdges ] = binVarByDist( inArc_NH.dist,inArc_NH.NC,binEdges );

    totesComboPlot( dist_inRow,NC_inRow,SEM_inRow,num_ROIs_inRow,NC_inRow,dist_inArc,NC_inArc,SEM_inArc,num_ROIs_inArc,NC_inArc )
    xlabel('distance')
ylabel('NC')
title('NH')
legend('in row','in arc')

%%

binEdges=10
[ dist_inRow,NC_inRow,SEM_inRow,num_ROIs_inRow,binEdges ] = binVarByDist( inRow_EN.dist,inRow_EN.NC,binEdges );
    [ dist_inArc,NC_inArc,SEM_inArc,num_ROIs_inArc,binEdges ] = binVarByDist( inArc_EN.dist,inArc_EN.NC,binEdges );

    totesComboPlot( dist_inRow,NC_inRow,SEM_inRow,num_ROIs_inRow,NC_inRow,dist_inArc,NC_inArc,SEM_inArc,num_ROIs_inArc,NC_inArc )
    xlabel('distance')
ylabel('NC')
title('EN')
legend('in row','in arc')

%%

figure; hold on
plot_scatterRLine(inArc_NH.dist,inArc_NH.NC)
xlabel('in arc distance, NH')
ylabel('NC')

figure; hold on
plot_scatterRLine(inArc_NH.dist,inArc_NH.SC)
xlabel('in arc distance, NH')
ylabel('SC')

figure; hold on
plot_scatterRLine(inArc_EN.dist,inArc_EN.NC)
xlabel('in arc distance, EN')
ylabel('NC')

figure; hold on
plot_scatterRLine(inArc_EN.dist,inArc_EN.SC)
xlabel('in arc distance, EN')
ylabel('SC')

%%

figure; hold on
plot_scatterRLine(inRow_NH.dist,inRow_NH.NC)
xlabel('in Row distance, NH')
ylabel('NC')

figure; hold on
plot_scatterRLine(inRow_NH.dist,inRow_NH.SC)
xlabel('in Row distance, NH')
ylabel('SC')

figure; hold on
plot_scatterRLine(inRow_EN.dist,inRow_EN.NC)
xlabel('in Row distance, EN')
ylabel('NC')

figure; hold on
plot_scatterRLine(inRow_EN.dist,inRow_EN.SC)
xlabel('in Row distance, EN')
ylabel('SC')

%%
figure; hold on
plot_4cdfs(inRow_EN.NC,inArc_EN.NC,inRow_NH.NC,inArc_NH.NC)
legend('in row,EN','in Arc,EN','in row, NH','in arc, EN')
xlabel('NC')

figure; hold on
plot_4cdfs(inRow_EN.dist,inArc_EN.dist,inRow_NH.dist,inArc_NH.dist)
legend('in row,EN','in Arc,EN','in row, NH','in arc, NH')
xlabel('inter-ROI distances')

figure; hold on
plot_4cdfs(inRow_EN.SC,inArc_EN.SC,inRow_NH.SC,inArc_NH.SC)
legend('in row,EN','in Arc,EN','in row, NH','in arc, NH')
xlabel('SC')

end