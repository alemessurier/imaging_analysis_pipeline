function compare_SCvNC_inBarr(paths_NH,paths_EN,npSub )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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


%% plot SC vs NC within barrels, by field

% EN
for K=1:length(results_EN)
    NCs=mean(results_EN(K).NC_inBarr,2);
    SCs=results_EN(K).SC_inBarr;
    [r2_EN(K),p_EN(K),slope_inBarr_EN(K)]=plot_scatterRLine(SCs,NCs);
    xlabel('EN, SC')
    ylabel('EN, NC')
end

%NH
for K=1:length(results_NH)
    NCs=mean(results_NH(K).NC_inBarr,2);
    SCs=results_NH(K).SC_inBarr;
    [r2_NH(K),p_NH(K),slope_inBarr_NH(K)]=plot_scatterRLine(SCs,NCs);
    xlabel('NH, SC')
    ylabel('NH, NC')
end

[ pval_inBarr] = notBoxPlot_ENvNH( slope_inBarr_EN,slope_inBarr_NH )
title(['pairs in barrels, p=',num2str(pval_inBarr)])
ylabel('SC vs NC slope by field')
%% plot SC vs NC across barrels, by field

% EN
for K=1:length(results_EN)
    NCs=mean(results_EN(K).NC_xBarr,2);
    SCs=results_EN(K).SC_xBarr;
    if ~isempty(NCs)
        [r2_EN(K),p_EN(K),slope_xBarr_EN(K)]=plot_scatterRLine(SCs,NCs);
        xlabel('EN, SC')
        ylabel('EN, NC')
    else
        r2_EN(K)=nan;
        p_EN(K)=nan;
        slope_xBarr_EN(K)=nan;
    end
end

%NH
for K=1:length(results_NH)
    NCs=mean(results_NH(K).NC_xBarr,2);
    SCs=results_NH(K).SC_xBarr;
    if ~isempty(NCs)
        [r2_NH(K),p_NH(K),slope_xBarr_NH(K)]=plot_scatterRLine(SCs,NCs);
        xlabel('NH, SC')
        ylabel('NH, NC')
    else
        r2_EN(K)=nan;
        p_EN(K)=nan;
        slope_xBarr_NH(K)=nan;
    end
end

[ pval_xBarr] = notBoxPlot_ENvNH( slope_xBarr_EN,slope_xBarr_NH );
title(['pairs across barrels, p=',num2str(pval_xBarr)])
ylabel('SC vs NC slope by field')
%% compare slope of SC v NC across and within barrels, by field

EN_slopes=[slope_inBarr_EN',slope_xBarr_EN'];
NH_slopes=[slope_inBarr_NH',slope_xBarr_NH'];

figure; hold on

for K=1:length(NH_slopes)
    plot(1:2,NH_slopes(K,:),'ko-','LineWidth',1)
end

for K=1:length(EN_slopes)
    plot(1:2,EN_slopes(K,:),'ro-','LineWidth',1)
end

mean_slope_EN=nanmean(EN_slopes);
mean_slope_NH=nanmean(NH_slopes);

plot(1:2,mean_slope_EN,'ro-','LineWidth',3)
plot(1:2,mean_slope_NH,'ko-','LineWidth',3)

tmp=gca;
tmp.XLim=[0.5 2.5];
tmp.XTick=1:2;
tmp.XTickLabel={'in barrel','across barrel'}
ylabel('SC vs NC slope')


%% whole population
SCs_EN_inBarr=cat(1,results_EN(:).SC_inBarr);
NCs_EN_inBarr=cat(1,results_EN(:).NC_inBarr);
NCs_EN_inBarr=mean(NCs_EN_inBarr,2);
plot_scatterRLine(SCs_EN_inBarr,NCs_EN_inBarr)

SCs_EN_xBarr=cat(1,results_EN(:).SC_xBarr);
NCs_EN_xBarr=cat(1,results_EN(:).NC_xBarr);
NCs_EN_xBarr=mean(NCs_EN_xBarr,2);
plot_scatterRLine(SCs_EN_xBarr,NCs_EN_xBarr)


%% whole population
SCs_NH_inBarr=cat(1,results_NH(:).SC_inBarr);
NCs_NH_inBarr=cat(1,results_NH(:).NC_inBarr);
NCs_NH_inBarr=mean(NCs_NH_inBarr,2);
plot_scatterRLine(SCs_NH_inBarr,NCs_NH_inBarr)

SCs_NH_xBarr=cat(1,results_NH(:).SC_xBarr);
NCs_NH_xBarr=cat(1,results_NH(:).NC_xBarr);
NCs_NH_xBarr=mean(NCs_NH_xBarr,2);
plot_scatterRLine(SCs_NH_xBarr,NCs_NH_xBarr)

end

