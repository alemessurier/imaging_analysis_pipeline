function pvals=somatotopy_plots_2017( paths_NH,paths_EN,npSub )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

 %% colorcode each imaging field by BW
for K=1:length(paths_EN)
     colorCodeROIs_byWhiskPref( paths_EN{K});
     title(paths_EN{K})
end
for K=1:length(paths_NH)
    colorCodeROIs_byWhiskPref( paths_NH{K});
    title(paths_NH{K})
end

%% plot percent of cells tuned to CW by radius (EBW)

[mean_dist_NH,mean_percMatch_NH,sem_percMatch_NH,percMatch_byBarr_NH,percMatch_NH,num_ROIs_NH,CWtunedInds_binned_NH,binEdges]=plot_CWtunedByRad(paths_NH,npSub,10);
[mean_dist_EN,mean_percMatch_EN,sem_percMatch_EN,percMatch_byBarr_EN, percMatch_EN,num_ROIs_EN,CWtunedInds_binned_EN]=plot_CWtunedByRad(paths_EN,npSub,binEdges);

% figure; hold on
% for i=1:length(mean_dist_NH)
%     pl(i)=notBoxPlot(percMatch_byBarr_NH{i},mean_dist_NH(i),[],'line')
%     pl(i).sd.Color='none';
%     pl(i).data.MarkerFaceColor='none';
%     pl(i).data.Color=[0 0 0];
%     pl(i).data.MarkerSize=4;
%     pl(i).mu.MarkerSize=12;
%     pl(i).mu.Color=[0 0 0];
%     pl(i).mu.MarkerFaceColor=[0 0 0];
% end
% 
% 
% for i=1:length(mean_dist_EN)
%     pr(i)=notBoxPlot(percMatch_byBarr_EN{i},mean_dist_EN(i),[],'line')
%     pr(i).sd.Color='none';
%     pr(i).data.MarkerFaceColor='none';
%     pr(i).data.Color=[1 0 0];
%     pr(i).data.MarkerSize=4;
%     pr(i).mu.MarkerSize=12;
%     pr(i).mu.Color=[1 0 0];
% end

for i=1:numel(CWtunedInds_binned_NH)
    P(i)=permutationTest(CWtunedInds_binned_NH{i},CWtunedInds_binned_EN{i},10000);
end

figure;
ax1(1)=axes('XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none');
hold on
plot(mean_dist_NH,num_ROIs_NH,'k:');
plot(mean_dist_EN,num_ROIs_EN,'r:');

ylabel('# ROIs')


position=ax1(1).Position;
ax1(2) = axes('Position',position,...
    'Color','none');
ax1(2).XLim=ax1(1).XLim;
hold on
ax1(1).XTickLabel=[];

plot(mean_dist_NH,percMatch_NH,'ko-');
plot(mean_dist_EN,percMatch_EN,'ro-');
xlabel('distance from CW center (um)')
ylabel('% tuned to CW (EBW)')

legend('NH','EN')
linkaxes(ax1,'x')


% figure;
% 
% ax2(1)=axes('XAxisLocation','bottom',...
%     'YAxisLocation','right',...
%     'Color','none');
% hold on
% plot(mean_dist_NH,num_ROIs_NH,'k:');
% plot(mean_dist_EN,num_ROIs_EN,'r:');
% 
% ylabel('# ROIs')
% 
% 
% position=ax2(1).Position;
% ax2(2) = axes('Position',position,...
%     'Color','none');
% ax2(2).XLim=ax2(1).XLim;
% hold on
% 
% 
% [hl,hp]=boundedline(mean_dist_NH,mean_percMatch_NH,sem_percMatch_NH,'ko-',mean_dist_EN,mean_percMatch_EN,sem_percMatch_EN,'ro-','alpha')
% xlabel('distance from CW center (um)')
% ylabel('mean % of cells tuned to CW (EBW)')
% hp_info=get(hp);
% hp_ingo(1).Annotation.LegendInformation.IconDisplayStyle='off';
% hp_ingo(2).Annotation.LegendInformation.IconDisplayStyle='off';
% legend('NH','EN')
% ax2(1).XTickLabel=[];
% linkaxes(ax2,'x')
% % for i=1:10
% %     [~,pval]=ttest2(percMatch_byBarr_NH{i},percMatch_byBarr_EN{i});
% %     %     if pval<0.05
% %     %         plot(mean_dist_NH(i),mean_percMatch_NH(i)+2*sem_percMatch_NH(i),'r*')
% %     %     end
% %     pvals(i)=pval;
% % end


%% plot percent of cells tuned to CW by radius (absolute BW)
[mean_dist_num_NH,mean_percMatch_num_NH,sem_percMatch_num_NH,percMatch_byBarr_num_NH,percMatch_num_NH,num_ROIs_num_NH,binEdges]=plot_CWtunedByRad_num(paths_NH,npSub,10);
[mean_dist_num_EN,mean_percMatch_num_EN,sem_percMatch_num_EN,percMatch_byBarr_num_EN, percMatch_num_EN,num_ROIs_num_EN]=plot_CWtunedByRad_num(paths_EN,npSub,binEdges);




figure;
ax1(1)=axes('XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none');
hold on
plot(mean_dist_num_NH,num_ROIs_num_NH,'k:');
plot(mean_dist_num_EN,num_ROIs_num_EN,'r:');

ylabel('# ROIs')


position=ax1(1).Position;
ax1(2) = axes('Position',position,...
    'Color','none');
ax1(2).XLim=ax1(1).XLim;
hold on
ax1(1).XTickLabel=[];

plot(mean_dist_num_NH,percMatch_num_NH,'ko-');
plot(mean_dist_num_EN,percMatch_num_EN,'ro-');
xlabel('distance from CW center (um)')
ylabel('% tuned to CW')

legend('NH','EN')
linkaxes(ax1,'x')


% figure;
% 
% ax2(1)=axes('XAxisLocation','bottom',...
%     'YAxisLocation','right',...
%     'Color','none');
% hold on
% plot(mean_dist_num_NH,num_ROIs_num_NH,'k:');
% plot(mean_dist_num_EN,num_ROIs_num_EN,'r:');
% 
% ylabel('# ROIs')
% 
% 
% position=ax2(1).Position;
% ax2(2) = axes('Position',position,...
%     'Color','none');
% ax2(2).XLim=ax2(1).XLim;
% hold on
% 
% 
% [hl,hp]=boundedline(mean_dist_num_NH,mean_percMatch_num_NH,sem_percMatch_num_NH,'ko-',mean_dist_num_EN,mean_percMatch_num_EN,sem_percMatch_num_EN,'ro-','alpha')
% xlabel('distance from CW center (um)')
% ylabel('mean % of cells tuned to CW')
% hp_info=get(hp);
% hp_ingo(1).Annotation.LegendInformation.IconDisplayStyle='off';
% hp_ingo(2).Annotation.LegendInformation.IconDisplayStyle='off';
% legend('NH','EN')
% ax2(1).XTickLabel=[];
% linkaxes(ax2,'x')
% 

%%
[mean_dist_NH,mean_responseZ_NH,sem_responseZ_NH,mean_responseDF_NH,sem_responseDF_NH,num_ROIs_NH,CW_Z_NH,CW_dF_NH,responsesZ_NH,responsesDF_NH,binEdges]=plot_CWresponseByRad(paths_NH,npSub,10);
[mean_dist_EN,mean_responseZ_EN,sem_responseZ_EN,mean_responseDF_EN,sem_responseDF_EN,num_ROIs_EN,CW_Z_EN,CW_dF_EN,responsesZ_EN,responsesDF_EN]=plot_CWresponseByRad(paths_EN,npSub,binEdges);

%% plots

% Zscored response by distance
indsENplot=num_ROIs_EN>0;
mean_dist_EN=mean_dist_EN(indsENplot);
mean_responseZ_EN=mean_responseZ_EN(indsENplot);
sem_responseZ_EN=sem_responseZ_EN(indsENplot);
num_ROIs_EN=num_ROIs_EN(indsENplot);

[Z_respByDist,ax]=totesComboPlot( mean_dist_NH,mean_responseZ_NH,...
    sem_responseZ_NH,num_ROIs_NH,CW_Z_NH,mean_dist_EN,mean_responseZ_EN,sem_responseZ_EN,num_ROIs_EN,CW_Z_EN )

ax(2).YLabel.String='Z-scored whisker response';
ax(2).XLabel.String='distance from whisker column';

for i=1:10
    pval=permutationTest_mean(responsesDF_NH{i},responsesDF_EN{i},10000);
%     if pval<0.05 && pval>0.01
%         plot(mean_dist_NH(i),0,'r*')
%     elseif pval<0.01 && pval>0.001
%         plot(mean_dist_NH(i),0,'r+')
%     elseif pval<0.001
%         plot(mean_dist_NH(i),0,'rsquare')
%     end
    
    pvals(i)=pval;
end

% median dF/F response by distance

mean_responseDF_EN=mean_responseDF_EN(indsENplot);
sem_responseDF_EN=sem_responseDF_EN(indsENplot);


[dF_respByDist,ax]=totesComboPlot( mean_dist_NH,mean_responseDF_NH,...
    sem_responseDF_NH,num_ROIs_NH,CW_dF_NH,mean_dist_EN,mean_responseDF_EN,sem_responseDF_EN,num_ROIs_EN,CW_dF_EN )

ax(2).YLabel.String='median dF/F';
ax(2).XLabel.String='distance from whisker column';



% %%
% for K=1:length(paths_EN)
%     switch npSub
%         case 0
%             load(strcat(paths_EN{K},'whiskPref.mat'),'whiskPref');
%             cd(paths_EN{K})
%             fname_s1=dir('step1_*');
%             load(strcat(paths_EN{K},fname_s1(end).name),'permTestResults','traceByStim');
%         case 1
%             load(strcat(paths_EN{K},'whiskPref_NPsub.mat'),'whiskPref');
%     end
%     cellNames=fieldnames(whiskPref);
%     numBW_EN{K}=cellfun(@(x)length(whiskPref.(x){1}),cellNames,'Uni',1);
%     [ ~,sig_inds_whisk ] = find_sigROIs( permTestResults,traceByStim );
%     numSigW_EN{K}=cellfun(@sum,sig_inds_whisk);
% end
% 
% numBW_EN_all=vertcat(numBW_EN{:});
% numSigW_EN_all=horzcat(numSigW_EN{:});
% 
% for K=1:length(paths_NH)
%     switch npSub
%         case 0
%             load(strcat(paths_NH{K},'whiskPref.mat'),'whiskPref');
%             cd(paths_NH{K})
%             fname_s1=dir('step1_*');
%             load(strcat(paths_NH{K},fname_s1(end).name),'permTestResults','traceByStim');
%         case 1
%             load(strcat(paths_NH{K},'whiskPref_NPsub.mat'),'whiskPref');
%     end
%     cellNames=fieldnames(whiskPref);
%     numBW_NH{K}=cellfun(@(x)length(whiskPref.(x){1}),cellNames,'Uni',1);
%     [ ~,sig_inds_whisk ] = find_sigROIs( permTestResults,traceByStim );
%     numSigW_NH{K}=cellfun(@sum,sig_inds_whisk);
% end
% numBW_NH_all=vertcat(numBW_NH{:});
% numSigW_NH_all=horzcat(numSigW_NH{:});
% 
% figure; hold on
% cdfplot(numBW_EN_all)
% cdfplot(numBW_NH_all)
% legend('EN','NH')
% ylabel('')
% xlabel('# equivalent BWs')
% [~,pval]=ttest2(numBW_EN_all,numBW_NH_all);
% title(strcat('p=',num2str(pval)))
% 
% % figure; hold on
% % hist(numBW_EN_all,1:9);
% % ax=gca;
% % ax.XTick=1:9;
% % ax.XTickLabel=1:9;
% % ax.Children.FaceColor='none';
% % ax.Children.LineWidth=1.5;
% % ax.Children.EdgeColor='b';
% % xlabel('# equivalent BWs')
% % ylabel('number of ROIs')
% % title('EN')
% 
% figure; hold on
% histogram(numBW_NH_all,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
% histogram(numBW_EN_all,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
% ax=gca;
% ax.XTick=1:9;
% ax.XTickLabel=1:9;
% title(strcat('p=',num2str(pval)))
% 
% % ax.Children.FaceColor='none';
% % ax.Children(1).LineWidth=1.5;
% % ax.Children(2).LineWidth=1.5;
% 
% xlabel('# equivalent BWs')
% ylabel('fraction of ROIs')
% 
% 
% figure; hold on
% histogram(numSigW_NH_all,1:11,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
% histogram(numSigW_EN_all,1:11,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
% ax=gca;
% ax.XTick=1:10;
% ax.XTickLabel=0:9;
% % ax.Children.FaceColor='none';
% % ax.Children(1).LineWidth=1.5;
% % ax.Children(2).LineWidth=1.5;
% %
% xlabel('# whiskers in RF')
% ylabel('fraction of ROIs')
% 
% [~,p]=ttest2(numSigW_NH_all,numSigW_EN_all)
% title(['p=',num2str(p)])
% 
% %%
% 
% [percentBWmatchCW_EN,percentBWmatchCW_num_EN]=plot_BWmatchCW(paths_EN);
% [percentBWmatchCW_NH,percentBWmatchCW_num_NH]=plot_BWmatchCW(paths_NH);
% 
% figure; hold on
% 
% pl(1)=notBoxPlot(percentBWmatchCW_EN,1,[],'line')
% pl(1).sd.Color='none';
% pl(1).data.MarkerFaceColor='none';
% 
% 
% pl(2)=notBoxPlot(percentBWmatchCW_NH,2,[],'line')
% pl(2).sd.Color='none';
% pl(2).data.MarkerFaceColor='none';
% set(gca,'XTick',[1:2])
% set(gca,'XTickLabels',{'EN','NH'})
% set(gca,'XLim',[0 3])
% set(gca,'FontSize',18)
% [~,p]=ttest2(percentBWmatchCW_EN,percentBWmatchCW_NH);
% ylabel('fraction of matched cells/barrel')
% % title(strcat('p=',num2str(p)))
% title(['CW is statistically equivalent to BW, p=',num2str(p)],'FontSize',18,'FontWeight','normal')
% 
% figure; hold on
% 
% pl(1)=notBoxPlot(percentBWmatchCW_num_EN,1,[],'line')
% pl(1).sd.Color='none';
% pl(1).data.MarkerFaceColor='none';
% 
% 
% pl(2)=notBoxPlot(percentBWmatchCW_num_NH,2,[],'line')
% pl(2).sd.Color='none';
% pl(2).data.MarkerFaceColor='none';
% set(gca,'XTick',[1:2])
% set(gca,'XTickLabels',{'EN','NH'})
% set(gca,'XLim',[0 3])
% set(gca,'FontSize',18)
% [~,p]=ttest2(percentBWmatchCW_num_EN,percentBWmatchCW_num_NH);
% ylabel('fraction of matched cells/barrel')
% % title(strcat('p=',num2str(p)))
% title(['CW=BW, p=',num2str(p)],'FontSize',18,'FontWeight','normal')

%%
    function [mean_dist,mean_responseZ,sem_responseZ,mean_responseDF,...
            sem_responseDF,num_ROIs,CW_response_Z,CW_response_dF,all_responseZ,all_responseDF,binEdges]=plot_CWresponseByRad(pathNames,npSub,binEdges)
        
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
                    traceByStim=traceByStim(2);
                    sponTrace=sponTrace(2);
                    permTestResults=permTestResults(2);
                    whiskPref=whiskPref(2);
            end
            sigROIs=find_sigROIs(permTestResults,traceByStim);
            [results  ] = find_CWresponse_byRad(sigROIs,traceByStim,sponTrace,framesEvoked,ROItoBarrel);
            results_all=[results_all,results];
            
        end
        
        
        
        distsAll=cat(1,results_all(:).ROI_dists);
        response_Z_all=cat(2,results_all(:).response_Z);
        response_dF_all=cat(2,results_all(:).response_dF);
        CW_response_Z=response_Z_all(distsAll<160);
        CW_response_dF=response_dF_all(distsAll<160);
        
        [distsAll,sortInds]=sort(distsAll,'ascend');
        response_Z_all=response_Z_all(sortInds);
        response_dF_all=response_dF_all(sortInds);
        [~,binEdges,binInds]=histcounts(distsAll,binEdges);
        
        
        for i=1:10
            all_responseZ{i}=response_Z_all(binInds==i);
            all_responseDF{i}=response_dF_all(binInds==i);
            mean_responseZ(i)=nanmean(response_Z_all(binInds==i));
            sem_responseZ(i)=nanstd(response_Z_all(binInds==i))/sqrt(sum(binInds==i));
            
            mean_responseDF(i)=nanmean(response_dF_all(binInds==i));
            sem_responseDF(i)=nanstd(response_dF_all(binInds==i))/sqrt(sum(binInds==i));
            
            mean_dist(i)=nanmean(distsAll(binInds==i));
            num_ROIs(i)=sum(binInds==i);
        end
    end


    function [mean_dist,mean_percMatch,sem_percMatch,percMatch_byBarr,percMatch,num_ROIs,binEdges]=plot_CWtunedByRad_num(pathNames,npSub,binEdges)
        
        
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
                    traceByStim=traceByStim(2);
                    sponTrace=sponTrace(2);
                    permTestResults=permTestResults(2);
                    whiskPref=whiskPref(2);
            end
            sigROIs=find_sigROIs(permTestResults,traceByStim);
            results=find_CWtuned_byRad_num(traceByStim,sigROIs,framesEvoked,ROItoBarrel);
            results_all=[results_all,results];
            
        end
        
        
        
        distsAll=cat(1,results_all(:).ROI_dists);
        percMatch_all=cat(2,results_all(:).percentMatch);
        
        CWtunedInds=cat(1,results_all(:).CWtunedInds);
        [distsAll,sortInds]=sort(distsAll,'ascend');
        CWtunedInds=CWtunedInds(sortInds);
        percMatch_all=percMatch_all(sortInds);
        
        [~,binEdges,binInds]=histcounts(distsAll,binEdges);
        
        for i=1:10
            mean_percMatch(i)=nanmean(percMatch_all(binInds==i));
            sem_percMatch(i)=nanstd(percMatch_all(binInds==i))/sqrt(sum(binInds==i));
            mean_dist(i)=nanmean(distsAll(binInds==i));
            percMatch_byBarr{i}=percMatch_all(binInds==i);
            percMatch(i)=sum(CWtunedInds(binInds==i))/sum(binInds==i);
            num_ROIs(i)=sum(binInds==i);
        end
    end

    function [mean_dist,mean_percMatch,sem_percMatch,percMatch_byBarr,percMatch,num_ROIs,CWtunedInds_binned,binEdges]=plot_CWtunedByRad(pathNames,npSub,binEdges)
        
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
                    traceByStim=traceByStim(2);
                    sponTrace=sponTrace(2);
                    permTestResults=permTestResults(2);
                    whiskPref=whiskPref(2);
            end
            sigROIs=find_sigROIs(permTestResults,traceByStim);
            results=find_CWtuned_byRad( traceByStim, sigROIs,whiskPref,ROItoBarrel);
            results_all=[results_all,results];
            
        end
        
        
        
        distsAll=cat(1,results_all(:).ROI_dists);
        percMatch_all=cat(2,results_all(:).percentMatch);
        
        CWtunedInds=cat(1,results_all(:).CWtunedInds);
        [distsAll,sortInds]=sort(distsAll,'ascend');
        CWtunedInds=CWtunedInds(sortInds);
        percMatch_all=percMatch_all(sortInds);
        [~,binEdges,binInds]=histcounts(distsAll,binEdges);
        
        
        for i=1:10
            mean_percMatch(i)=nanmean(percMatch_all(binInds==i));
            sem_percMatch(i)=nanstd(percMatch_all(binInds==i))/sqrt(sum(binInds==i));
            mean_dist(i)=nanmean(distsAll(binInds==i));
            percMatch_byBarr{i}=percMatch_all(binInds==i);
            percMatch(i)=sum(CWtunedInds(binInds==i))/sum(binInds==i);
            CWtunedInds_binned{i}=CWtunedInds(binInds==i);
            num_ROIs(i)=sum(binInds==i);
        end
    end

    function [percentBWmatchCW,percentBWmatchCW_num]=plot_BWmatchCW(pathNames)
        
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
                    traceByStim=traceByStim(2);
                    sponTrace=sponTrace(2);
                    permTestResults=permTestResults(2);
                    whiskPref=whiskPref(2);
            end
            [percentBWmatch{J},percentBWmatch_num{J}] = find_BWmatchCW( ROIsInBarrel,whiskPref,traceByStim,framesEvoked );
        end
        percentBWmatchCW=horzcat(percentBWmatch{:});
        percentBWmatchCW_num=horzcat(percentBWmatch_num{:});
        
    end
end