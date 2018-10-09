function  [fullfield_NH,fullfield_EN,inBarr_NH,inBarr_EN,xBarr_NH,xBarr_EN]=make_CorrAnalysisPlots( paths_NH,paths_EN,npSub,rval,ROIs_NH,ROIs_EN )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% find NC,SC for all ROIs in field - enriched
results_EN=[];
for K=1:length(paths_EN)
    switch npSub
        case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_EN{K} );
        case 1
            [traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L23( paths_EN{K},rval );
            
    end
    
    if isempty(ROIs_EN)
        sigROIs=find_sigROIs(permTestResults,traceByStim);
    else
        sigROIs=ROIs_EN{K};
    end
    [ results ] = find_corr_sameBW( traceByStim,sponTrace,sigROIs,framesEvoked,ROI_positions,mag);
    if ~isempty(results(1).NC)
        results_EN=[results_EN,results];
    end
end

%% find NC,SC for all ROIs in field - NH
results_NH=[];
for K=1:length(paths_NH)
    switch npSub
        case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_NH{K} );
        case 1
            [traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L23( paths_NH{K},rval );
            
    end
    
    
    if isempty(ROIs_NH)
        sigROIs=find_sigROIs(permTestResults,traceByStim);
    else
        sigROIs=ROIs_NH{K};
    end
    
    [ results ] = find_corr_sameBW( traceByStim,sponTrace,sigROIs,framesEvoked,ROI_positions,mag);
    if ~isempty(results(1).NC)
        results_NH=[results_NH,results];
    end
    
end


%% make plots of NC and SC by distance, whole field

NC_NH_all=cat(1,results_NH(:).NC);
fullfield_NH.NC=mean(NC_NH_all(:,1:9),2);
shuffNC_NH_all=cat(1,results_NH(:).shuffNC);
fullfield_NH.NC_shuff=mean(shuffNC_NH_all(:,1:9),2);

fullfield_NH.dist=cat(2,results_NH(:).ROIdistance);
fullfield_NH.SC=cat(1,results_NH(:).SC);
fullfield_NH.SC_shuff=cat(2,results_NH(:).shuffSC)';




%%


NC_EN_all=cat(1,results_EN(:).NC);
fullfield_EN.NC=mean(NC_EN_all(:,1:9),2);
shuffNC_EN_all=cat(1,results_EN(:).shuffNC);
fullfield_EN.NC_shuff=mean(shuffNC_EN_all(:,1:9),2);

fullfield_EN.dist=cat(2,results_EN(:).ROIdistance);
fullfield_EN.SC=cat(1,results_EN(:).SC);
fullfield_EN.SC_shuff=cat(2,results_EN(:).shuffSC)';


%
% for i=1:length(allNC_NH)
%     tmp_p(i)=permutationTest_mean(allNC_NH{i},allNC_EN{i},10000);
%     tmp_p2(i)=permutationTest_mean(allSC_NH{i},allSC_EN{i},10000);
% end
% pvals_NC.field=tmp_p;
% pvals_SC.field=tmp_p2;


%% make plots comparing SC,NC by distance within column pairs vs all

% find NC,SC for ROIs within and across columns
results_NH=[];
for K=1:length(paths_NH)
    switch npSub
        case 0
            [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_nonNPsub_data( paths_NH{K} );
        case 1
            [traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L23( paths_NH{K},rval );
            
    end
    
    
    if isempty(ROIs_NH)
            sigROIs=find_sigROIs(permTestResults,traceByStim);
    else
            sigROIs=ROIs_NH{K};
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
            [traceByStim,sponTrace,framesEvoked,permTestResults,...
                dists,ROIsInBarrel,ROItoBarrel,ROI_positions,~,~,~,~,mag ] = load_NPsub_data_L23( paths_EN{K},rval );
            
    end
    
    
    if isempty(ROIs_EN)
            sigROIs=find_sigROIs(permTestResults,traceByStim);
    else
            sigROIs=ROIs_EN{K};
    end
    
    [ results ] = find_corrInBarrelvsAcross( sigROIs,traceByStim,framesEvoked,ROIsInBarrel,sponTrace,ROI_positions,mag);
    if ~isempty(results)
        results_EN=[results_EN,results];
    end
end





%%
NC_EN_inBarr=cat(1,results_EN(:).NC_inBarr);
inBarr_EN.NC=mean(NC_EN_inBarr(:,1:9),2);
NC_EN_xBarr=cat(1,results_EN(:).NC_xBarr);
xBarr_EN.NC=mean(NC_EN_xBarr(:,1:9),2);

inBarr_EN.SC=cat(1,results_EN(:).SC_inBarr);
xBarr_EN.SC=cat(1,results_EN(:).SC_xBarr);
inBarr_EN.dist=cat(2,results_EN(:).ROIdist_inBarr);
xBarr_EN.dist=cat(2,results_EN(:).ROIdist_xBarr);

NC_NH_inBarr=cat(1,results_NH(:).NC_inBarr);
inBarr_NH.NC=mean(NC_NH_inBarr(:,1:9),2);
NC_NH_xBarr=cat(1,results_NH(:).NC_xBarr);
xBarr_NH.NC=mean(NC_NH_xBarr(:,1:9),2);

inBarr_NH.SC=cat(1,results_NH(:).SC_inBarr);
xBarr_NH.SC=cat(1,results_NH(:).SC_xBarr);
inBarr_NH.dist=cat(2,results_NH(:).ROIdist_inBarr);
xBarr_NH.dist=cat(2,results_NH(:).ROIdist_xBarr);

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

