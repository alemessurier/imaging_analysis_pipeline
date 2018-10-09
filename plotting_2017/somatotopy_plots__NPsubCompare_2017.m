function somatotopy_plots__NPsubCompare_2017( paths_NH,paths_EN)

%  %% colorcode each imaging field by BW
% for K=1:length(paths_EN)
%      colorCodeROIs_byWhiskPref_L4( paths_EN{K});
%      title(paths_EN{K})
% end
% for K=1:length(paths_NH)
% colorCodeROIs_byWhiskPref_L4( paths_NH{K});
%      title(paths_NH{K})
% end



%% plot percent of cells tuned to CW by radius (EBW)
for r=2
    [mean_dist_NH,mean_percMatch_NH,sem_percMatch_NH,percMatch_byBarr_NH,percMatch_NH,num_ROIs_NH]=plot_CWtunedByRad(paths_NH,r);
    [mean_dist_EN,mean_percMatch_EN,sem_percMatch_EN,percMatch_byBarr_EN, percMatch_EN,num_ROIs_EN]=plot_CWtunedByRad(paths_EN,r);
    
%     save(['E:\Data\reduced\CWtuned_EBW_L4_r',num2str(r),'.mat'],...
%         'mean_dist_NH','percMatch_NH','num_ROIs_NH','mean_dist_EN',...
%         'percMatch_EN','num_ROIs_EN')
    
    figure; hold on
    for i=1:length(mean_dist_NH)
        pl(i)=notBoxPlot(percMatch_byBarr_NH{i},mean_dist_NH(i),[],'line')
        pl(i).sd.Color='none';
        pl(i).data.MarkerFaceColor='none';
        pl(i).data.Color=[0 0 0];
        pl(i).data.MarkerSize=4;
        pl(i).mu.MarkerSize=12;
        pl(i).mu.Color=[0 0 0];
        pl(i).mu.MarkerFaceColor=[0 0 0];
    end
    
    
    for i=1:length(mean_dist_EN)
        pr(i)=notBoxPlot(percMatch_byBarr_EN{i},mean_dist_EN(i),[],'line')
        pr(i).sd.Color='none';
        pr(i).data.MarkerFaceColor='none';
        pr(i).data.Color=[1 0 0];
        pr(i).data.MarkerSize=4;
        pr(i).mu.MarkerSize=12;
        pr(i).mu.Color=[1 0 0];
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
    
    
    figure;
    
    ax2(1)=axes('XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none');
    hold on
    plot(mean_dist_NH,num_ROIs_NH,'k:');
    plot(mean_dist_EN,num_ROIs_EN,'r:');
    
    ylabel('# ROIs')
    
    
    position=ax2(1).Position;
    ax2(2) = axes('Position',position,...
        'Color','none');
    ax2(2).XLim=ax2(1).XLim;
    hold on
    
    
    [hl,hp]=boundedline(mean_dist_NH,mean_percMatch_NH,sem_percMatch_NH,'ko-',mean_dist_EN,mean_percMatch_EN,sem_percMatch_EN,'ro-','alpha')
    xlabel('distance from CW center (um)')
    ylabel('mean % of cells tuned to CW (EBW)')
    hp_info=get(hp);
    hp_ingo(1).Annotation.LegendInformation.IconDisplayStyle='off';
    hp_ingo(2).Annotation.LegendInformation.IconDisplayStyle='off';
    legend('NH','EN')
    ax2(1).XTickLabel=[];
    linkaxes(ax2,'x')
end
%% plot percent of cells tuned to CW by radius (absolute BW)
for r=2
    [mean_dist_num_NH,mean_percMatch_num_NH,sem_percMatch_num_NH,percMatch_byBarr_num_NH,percMatch_num_NH,num_ROIs_num_NH]=plot_CWtunedByRad_num(paths_NH,r);
    [mean_dist_num_EN,mean_percMatch_num_EN,sem_percMatch_num_EN,percMatch_byBarr_num_EN, percMatch_num_EN,num_ROIs_num_EN]=plot_CWtunedByRad_num(paths_EN,r);
    
%     save(['E:\Data\reduced\CWtuned_num_L4_r',num2str(r),'.mat'],...
%         'mean_dist_num_NH','percMatch_num_NH','num_ROIs_num_NH','mean_dist_num_EN',...
%         'percMatch_num_EN','num_ROIs_num_EN')
    
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
    title(num2str(r));
    
    figure;
    
    ax2(1)=axes('XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none');
    hold on
    plot(mean_dist_num_NH,num_ROIs_num_NH,'k:');
    plot(mean_dist_num_EN,num_ROIs_num_EN,'r:');
    
    ylabel('# ROIs')
    
    
    position=ax2(1).Position;
    ax2(2) = axes('Position',position,...
        'Color','none');
    ax2(2).XLim=ax2(1).XLim;
    hold on
    
    
    [hl,hp]=boundedline(mean_dist_num_NH,mean_percMatch_num_NH,sem_percMatch_num_NH,'ko-',mean_dist_num_EN,mean_percMatch_num_EN,sem_percMatch_num_EN,'ro-','alpha')
    xlabel('distance from CW center (um)')
    ylabel('mean % of cells tuned to CW')
    hp_info=get(hp);
    hp_ingo(1).Annotation.LegendInformation.IconDisplayStyle='off';
    hp_ingo(2).Annotation.LegendInformation.IconDisplayStyle='off';
    legend('NH','EN')
    ax2(1).XTickLabel=[];
    linkaxes(ax2,'x')
    title(num2str(r));
end

%%
for r=2
    [mean_dist_NH,mean_responseZ_NH,sem_responseZ_NH,mean_responseDF_NH,sem_responseDF_NH,num_ROIs_NH,CW_Z_NH,CW_dF_NH]=plot_CWresponseByRad(paths_NH,r);
    [mean_dist_EN,mean_responseZ_EN,sem_responseZ_EN,mean_responseDF_EN,sem_responseDF_EN,num_ROIs_EN,CW_Z_EN,CW_dF_EN]=plot_CWresponseByRad(paths_EN,r);
    %% plots
    
    % Zscored response by distance
    
    [Z_respByDist,ax]=totesComboPlot( mean_dist_NH,mean_responseZ_NH,...
        sem_responseZ_NH,num_ROIs_NH,CW_Z_NH,mean_dist_EN,mean_responseZ_EN,sem_responseZ_EN,num_ROIs_EN,CW_Z_EN )
    
    ax(2).YLabel.String='Z-scored whisker response';
    ax(2).XLabel.String='distance from whisker column';
    title(num2str(r));
    % median dF/F response by distance
    
    [dF_respByDist,ax]=totesComboPlot( mean_dist_NH,mean_responseDF_NH,...
        sem_responseDF_NH,num_ROIs_NH,CW_dF_NH,mean_dist_EN,mean_responseDF_EN,sem_responseDF_EN,num_ROIs_EN,CW_dF_EN )
    
    ax(2).YLabel.String='median dF/F';
    ax(2).XLabel.String='distance from whisker column';
    title(num2str(r));
end
%%
% %
% for r=2
%     numSigW_NH={};
%     numSigW_EN={};
%     for K=1:length(paths_EN)
%         cd(paths_EN{K})
%         fname_s1=dir('step1NPcorr_0*');
%         load(strcat(paths_EN{K},fname_s1(end).name),'permTestResults','traceByStim');
%         
%         [ ~,sig_inds_whisk ] = find_sigROIs( permTestResults(r),traceByStim(r) );
%         numSigW_EN{K}=cellfun(@sum,sig_inds_whisk);
%     end
%     
%     numSigW_EN_all=horzcat(numSigW_EN{:});
%     
%     for K=1:length(paths_NH)
%         cd(paths_NH{K})
%         fname_s1=dir('step1NPcorr_0*');
%         load(strcat(paths_NH{K},fname_s1(end).name),'permTestResults','traceByStim');
%         
%         [ ~,sig_inds_whisk ] = find_sigROIs( permTestResults(r),traceByStim(r) );
%         numSigW_NH{K}=cellfun(@sum,sig_inds_whisk);
%     end
%     
%     numSigW_NH_all=horzcat(numSigW_NH{:});
%     
%     
%     
%     figure; hold on
%     histogram(numSigW_NH_all,1:11,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
%     histogram(numSigW_EN_all,1:11,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
%     ax=gca;
%     ax.XTick=1:10;
%     ax.XTickLabel=0:9;
%     ax.Children.FaceColor='none';
%     ax.Children(1).LineWidth=1.5;
%     ax.Children(2).LineWidth=1.5;
%     
%     xlabel('# whiskers in RF')
%     ylabel('fraction of ROIs')
%     
%     [~,p]=ttest2(numSigW_NH_all,numSigW_EN_all);
%     title(['p=',num2str(p),', r=',num2str(r)])
%     
%     numSigW_NH_resp=numSigW_NH_all(numSigW_NH_all>0);
%     numSigW_EN_resp=numSigW_EN_all(numSigW_EN_all>0);
%     
%     figure; hold on
%     histogram(numSigW_NH_resp,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
%     histogram(numSigW_EN_resp,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
%     ax=gca;
%     ax.XTick=1:9;
%     ax.XTickLabel=1:9;
%     ax.Children.FaceColor='none';
%     ax.Children(1).LineWidth=1.5;
%     ax.Children(2).LineWidth=1.5;
%     
%     xlabel('# whiskers in RF')
%     ylabel('fraction of ROIs')
%     
%     [~,p]=ttest2(numSigW_NH_resp,numSigW_EN_resp);
%     title(['p=',num2str(p),', r=',num2str(r)])
% end

%%
    function [mean_dist,mean_responseZ,sem_responseZ,mean_responseDF,...
            sem_responseDF,num_ROIs,CW_response_Z,CW_response_dF]=plot_CWresponseByRad(pathNames,rval)
        
        results_all=[];
        for J=1:length(pathNames)
            
            [~,traceByStim,sponTrace,framesEvoked,permTestResults,...
                ~,~,ROItoBarrel] = load_NPsub_data( pathNames{J} );
            
            
            permTestResults=permTestResults(rval);
            traceByStim=traceByStim(rval);
            sponTrace=sponTrace(rval);
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
        [~,binEdges,binInds]=histcounts(distsAll,10);
        
        
        for i=1:10
            mean_responseZ(i)=nanmean(response_Z_all(binInds==i));
            sem_responseZ(i)=nanstd(response_Z_all(binInds==i))/sqrt(sum(binInds==i));
            
            mean_responseDF(i)=nanmean(response_dF_all(binInds==i));
            sem_responseDF(i)=nanstd(response_dF_all(binInds==i))/sqrt(sum(binInds==i));
            
            mean_dist(i)=nanmean(distsAll(binInds==i));
            num_ROIs(i)=sum(binInds==i);
        end
    end

    function [mean_dist,mean_percMatch,sem_percMatch,percMatch_byBarr,percMatch,num_ROIs]=plot_CWtunedByRad(pathNames,rval)
        
        results_all=[];
        for J=1:length(pathNames)
            [~,traceByStim,~,framesEvoked,permTestResults,...
                ~,~,ROItoBarrel,~,~,~,~,whiskPref ] = load_NPsub_data( pathNames{J} );
            
            permTestResults=permTestResults(rval);
            traceByStim=traceByStim(rval);
            whiskPref=whiskPref(rval);
            results=find_CWtuned_byRad( traceByStim, permTestResults,whiskPref,ROItoBarrel);
            results_all=[results_all,results];
            
        end
        
        
        
        distsAll=cat(1,results_all(:).ROI_dists);
        percMatch_all=cat(2,results_all(:).percentMatch);
        
        CWtunedInds=cat(1,results_all(:).CWtunedInds);
        [distsAll,sortInds]=sort(distsAll,'ascend');
        CWtunedInds=CWtunedInds(sortInds);
        percMatch_all=percMatch_all(sortInds);
        [~,binEdges,binInds]=histcounts(distsAll,10);
        
        
        for i=1:10
            mean_percMatch(i)=nanmean(percMatch_all(binInds==i));
            sem_percMatch(i)=nanstd(percMatch_all(binInds==i))/sqrt(sum(binInds==i));
            mean_dist(i)=nanmean(distsAll(binInds==i));
            percMatch_byBarr{i}=percMatch_all(binInds==i);
            percMatch(i)=sum(CWtunedInds(binInds==i))/sum(binInds==i);
            num_ROIs(i)=sum(binInds==i);
        end
    end

    function [mean_dist,mean_percMatch,sem_percMatch,percMatch_byBarr,percMatch,num_ROIs]=plot_CWtunedByRad_num(pathNames,rval)
        
        results_all=[];
        for J=1:length(pathNames)
            
            [~,traceByStim,~,framesEvoked,permTestResults,...
                ~,~,ROItoBarrel ] = load_NPsub_data( pathNames{J} );
            
            permTestResults=permTestResults(rval);
            traceByStim=traceByStim(rval);
            results=find_CWtuned_byRad_num(traceByStim,permTestResults,framesEvoked,ROItoBarrel);
            results_all=[results_all,results];
            
        end
        
        
        
        distsAll=cat(1,results_all(:).ROI_dists);
        percMatch_all=cat(2,results_all(:).percentMatch);
        
        CWtunedInds=cat(1,results_all(:).CWtunedInds);
        [distsAll,sortInds]=sort(distsAll,'ascend');
        CWtunedInds=CWtunedInds(sortInds);
        percMatch_all=percMatch_all(sortInds);
        [~,binEdges,binInds]=histcounts(distsAll,10);
        
        
        for i=1:10
            mean_percMatch(i)=nanmean(percMatch_all(binInds==i));
            sem_percMatch(i)=nanstd(percMatch_all(binInds==i))/sqrt(sum(binInds==i));
            mean_dist(i)=nanmean(distsAll(binInds==i));
            percMatch_byBarr{i}=percMatch_all(binInds==i);
            percMatch(i)=sum(CWtunedInds(binInds==i))/sum(binInds==i);
            num_ROIs(i)=sum(binInds==i);
        end
    end




end