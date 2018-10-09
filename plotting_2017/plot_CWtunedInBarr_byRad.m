function plot_CWtunedInBarr_byRad(paths_EN,paths_NH,npSub)

[mean_dist_EN,percMatch_EN,num_ROIs_EN,percentCWtuned_EN,numROIsBin_EN]=find_CWtuned_inBarr(paths_EN,npSub);
[mean_dist_NH,percMatch_NH,num_ROIs_NH,percentCWtuned_NH,numROIsBin_NH]=find_CWtuned_inBarr(paths_NH,npSub);


figure; hold on
for K=1:numel(mean_dist_EN)
    plot(mean_dist_EN{K},percMatch_EN{K},'ro');%,'LineWidth',0.5)
end

for K=1:numel(mean_dist_NH)
    plot(mean_dist_NH{K},percMatch_NH{K},'ko')%,'LineWidth',0.5)
end


figure; hold on

pl(1)=notBoxPlot(percentCWtuned_EN,1,[],'line')
pl(1).sd.Color='none';
%pl(1).data.MarkerFaceColor='none';
pl(1).data.Marker='none';
numROIsBin_EN=numROIsBin_EN(~isnan(percentCWtuned_EN));
XData=pl(1).data.XData;
YData=pl(1).data.YData;
hold on
for i=1:length(XData)
    plot(XData(i),YData(i),'ko','MarkerSize',numROIsBin_EN(i)/5)
end

pl(2)=notBoxPlot(percentCWtuned_NH,2,[],'line')
pl(2).sd.Color='none';
pl(2).data.MarkerFaceColor='none';
pl(2).data.Marker='none';
numROIsBin_NH=numROIsBin_NH(~isnan(percentCWtuned_NH));
XData=pl(2).data.XData;
YData=pl(2).data.YData;
hold on
for i=1:length(XData)
    plot(XData(i),YData(i),'ko','MarkerSize',numROIsBin_NH(i)/5)
end


set(gca,'XTick',[1:2])
set(gca,'XTickLabels',{'EN','NH'})
set(gca,'XLim',[0 3])
set(gca,'FontSize',18)
[~,p]=ttest2(percentCWtuned_EN,percentCWtuned_NH);
ylabel('fraction of matched cells/barrel')
% title(strcat('p=',num2str(p)))
title(['CW=BW, p=',num2str(p)],'FontSize',18,'FontWeight','normal')

    function [mean_dist,percMatch,num_ROIs,percentCWtuned,num_ROIs_inc]=find_CWtuned_inBarr(pathNames,npSub)
        binEdges=[0:75:1350];
        results_all=[];
        for J=1:length(pathNames)
            switch npSub
                case 0
                    [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                        dists,ROIsInBarrel,ROItoBarrel,ROI_positions,Stimuli,...
                        deltaF,sampRate,whiskPref ] = load_nonNPsub_data( pathNames{J} );
                case 1
                    [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                        dists,ROIsInBarrel,ROItoBarrel,ROI_positions,Stimuli,...
                        deltaF,sampRate,whiskPref ] = load_NPsub_data( pathNames{J} );
            end
            results=find_CWtuned_byRad( traceByStim, permTestResults,whiskPref,ROItoBarrel);
             results_all=[results_all,results];
            
            distsAll=cat(1,results(:).ROI_dists);
            percMatch_all=cat(2,results(:).percentMatch);
            
            CWtunedInds=cat(1,results(:).CWtunedInds);
            [distsAll,sortInds]=sort(distsAll,'ascend');
            CWtunedInds=CWtunedInds(sortInds);
            percMatch_all=percMatch_all(sortInds);
            [~,~,binInds]=histcounts(distsAll,binEdges);
            binIDs=unique(binInds);
            
            for i=1:length(binIDs)
                mean_dist{J}(i)=nanmean(distsAll(binInds==i));
                percMatch{J}(i)=sum(CWtunedInds(binInds==i))/sum(binInds==i);
                num_ROIs{J}(i)=sum(binInds==i);
            end
            
        end
        
        
        
        for i=1:numel(results_all)
            inds_inbarr=logical(results_all(i).ROI_dists<200);
            if sum(inds_inbarr)<10
                percentCWtuned(i)=nan;
            else
                CWtunedInBarr=results_all(i).CWtunedInds(inds_inbarr);
                percentCWtuned(i)=sum(CWtunedInBarr)/length(CWtunedInBarr);
                num_ROIs_inc(i)=length(CWtunedInBarr);
            end
        end
        
    end
end