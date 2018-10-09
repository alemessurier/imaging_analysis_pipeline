function pvals=somatotopy_plots_2017_L4( paths_NH,paths_EN,npSub )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%  %% colorcode each imaging field by BW
% for K=1:length(paths_EN)
%      colorCodeROIs_byWhiskPref_L4( paths_EN{K});
%      title(paths_EN{K})
% end
% for K=1:length(paths_NH)
%     colorCodeROIs_byWhiskPref_L4( paths_NH{K});
%     title(paths_NH{K})
% end

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

pvals.CWtunedByDist_EBW=P;

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
[mean_dist_num_NH,mean_percMatch_num_NH,sem_percMatch_num_NH,percMatch_byBarr_num_NH,percMatch_num_NH,num_ROIs_num_NH,CWtunedInds_binned_NH,binEdges]=plot_CWtunedByRad_num(paths_NH,npSub,10);
[mean_dist_num_EN,mean_percMatch_num_EN,sem_percMatch_num_EN,percMatch_byBarr_num_EN, percMatch_num_EN,num_ROIs_num_EN,CWtunedInds_binned_EN]=plot_CWtunedByRad_num(paths_EN,npSub,binEdges);

for i=1:numel(CWtunedInds_binned_NH)
    P(i)=permutationTest(CWtunedInds_binned_NH{i},CWtunedInds_binned_EN{i},10000);
end

pvals.CWtunedByDist_num=P;

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

% compute anova

% for i=1:length(CWtunedInds_binned_EN)
%     distInds_EN{i}=i*ones(size(CWtunedInds_binned_EN{i}));
%     distInds_NH{i}=i*ones(size(CWtunedInds_binned_NH{i}));
% end
% 
% testV=[cat(1,CWtunedInds_binned_EN{:}); cat(1,CWtunedInds_binned_NH{:})];
% distInds_EN=cat(1,distInds_EN{:});
% distInds_NH=cat(1,distInds_NH{:});
% EN=cell(size(distInds_EN));
% EN(:)={'EN'};
% NH=cell(size(distInds_NH));
% NH(:)={'NH'};
% distInds=[distInds_EN; distInds_NH];
% groupID=[EN(:); NH(:)]
% 
% 
% [P,T,STATS,TERMS]=anovan(testV,{distInds,groupID},'varnames',{'distance','group'},'model','interaction')
% pvals.CWtunedByDist=P(3);


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
    pvalZ=permutationTest_mean(responsesZ_NH{i},responsesZ_EN{i},10000);
    %     if pval<0.05 && pval>0.01
    %         plot(mean_dist_NH(i),0,'r*')
    %     elseif pval<0.01 && pval>0.001
    %         plot(mean_dist_NH(i),0,'r+')
    %     elseif pval<0.001
    %         plot(mean_dist_NH(i),0,'rsquare')
    %     end
    
    P(i)=pval;
    PZ(i)=pvalZ;
end

% pvals.RbyDistDF=P;
% pvals.RbyDistZ=PZ;
% median dF/F response by distance

mean_responseDF_EN=mean_responseDF_EN(indsENplot);
sem_responseDF_EN=sem_responseDF_EN(indsENplot);


[dF_respByDist,ax]=totesComboPlot( mean_dist_NH,mean_responseDF_NH,...
    sem_responseDF_NH,num_ROIs_NH,CW_dF_NH,mean_dist_EN,mean_responseDF_EN,sem_responseDF_EN,num_ROIs_EN,CW_dF_EN )

ax(2).YLabel.String='median dF/F';
ax(2).XLabel.String='distance from whisker column';

for i=1:length(responsesDF_EN)
    distInds_EN{i}=i*ones(size(responsesDF_EN{i}));
    distInds_NH{i}=i*ones(size(responsesDF_NH{i}));
end

testV=[responsesDF_EN{:}, responsesDF_NH{:}];
distInds_EN=[distInds_EN{:}];
distInds_NH=[distInds_NH{:}];
EN=cell(size(distInds_EN));
EN(:)={'EN'};
NH=cell(size(distInds_NH));
NH(:)={'NH'};
distInds=[distInds_EN, distInds_NH];
groupID=[EN(:); NH(:)];

figure;
[P,T,STATS,TERMS]=anovan(testV,{distInds,groupID},'varnames',{'distance','group'},'model','interaction')
pvals.RbyDistDF=P(2);


testV=[responsesZ_EN{:}, responsesZ_NH{:}];
figure;
[P,T,STATS,TERMS]=anovan(testV,{distInds,groupID},'varnames',{'distance','group'},'model','interaction')
pvals.RbyDistZ=P(2);

%%
for K=1:length(paths_EN)
    switch npSub
        case 0
            [~,traceByStim,~,~,permTestResults,...
                ~,~,~,~,~,~,~,whiskPref ] = load_nonNPsub_data( paths_EN{K});
        case 1
            [~,traceByStim,~,~,permTestResults,...
                ~,~,~,~,~,~,~,whiskPref ] = load_NPsub_data_L4( paths_EN{K},2 );
    end
    [ cellNames,sig_inds_whisk ] = find_sigROIs( permTestResults,traceByStim );
    numBW_EN{K}=cellfun(@(x)length(whiskPref.(x){1}),cellNames,'Uni',1);
    
    numSigW_EN{K}=cellfun(@sum,sig_inds_whisk);
end

numBW_EN_all=vertcat(numBW_EN{:});
numSigW_EN_all=horzcat(numSigW_EN{:});

for K=1:length(paths_NH)
    switch npSub
        case 0
            [~,traceByStim,~,~,permTestResults,...
                ~,~,~,~,~,~,~,whiskPref ] = load_nonNPsub_data( paths_NH{K} );
        case 1
            [~,traceByStim,~,~,permTestResults,...
                ~,~,~,~,~,~,~,whiskPref ] = load_NPsub_data_L4( paths_NH{K},2 );
    end
    [ cellNames,sig_inds_whisk ] = find_sigROIs( permTestResults,traceByStim );
    numBW_NH{K}=cellfun(@(x)length(whiskPref.(x){1}),cellNames,'Uni',1);
   
    numSigW_NH{K}=cellfun(@sum,sig_inds_whisk);
end
numBW_NH_all=vertcat(numBW_NH{:});
numSigW_NH_all=horzcat(numSigW_NH{:});

figure; hold on
cdfplot(numBW_EN_all)
cdfplot(numBW_NH_all)
legend('EN','NH')
ylabel('')
xlabel('# equivalent BWs')
[~,pval]=ttest2(numBW_EN_all,numBW_NH_all);
title(strcat('p=',num2str(pval)))

% figure; hold on
% hist(numBW_EN_all,1:9);
% ax=gca;
% ax.XTick=1:9;
% ax.XTickLabel=1:9;
% ax.Children.FaceColor='none';
% ax.Children.LineWidth=1.5;
% ax.Children.EdgeColor='b';
% xlabel('# equivalent BWs')
% ylabel('number of ROIs')
% title('EN')

figure; hold on
histogram(numBW_NH_all,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(numBW_EN_all,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
ax=gca;
ax.XTick=1:9;
ax.XTickLabel=1:9;
title(strcat('p=',num2str(pval)))

% ax.Children.FaceColor='none';
% ax.Children(1).LineWidth=1.5;
% ax.Children(2).LineWidth=1.5;

xlabel('# equivalent BWs')
ylabel('fraction of ROIs')


figure; hold on
histogram(numSigW_NH_all,1:11,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(numSigW_EN_all,1:11,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
ax=gca;
ax.XTick=1:10;
ax.XTickLabel=0:9;
% ax.Children.FaceColor='none';
% ax.Children(1).LineWidth=1.5;
% ax.Children(2).LineWidth=1.5;
%
xlabel('# whiskers in RF')
ylabel('fraction of ROIs')

[~,p]=ttest2(numSigW_NH_all,numSigW_EN_all)
title(['p=',num2str(p)])

%% compare RF sharpness

[rankedRFs_NH,rankedRFs_CWnorm_NH,rankedRFs_sponNorm_NH]=plot_rankedRFs(paths_NH,npSub)
[rankedRFs_EN,rankedRFs_CWnorm_EN,rankedRFs_sponNorm_EN]=plot_rankedRFs(paths_EN,npSub)

f1=figure; hold on
f2=figure; hold on
f3=figure; hold on
for K=1:length(paths_NH)
    figure(f1)
    plot(1:9,mean(rankedRFs_NH{K},1),'k.-')
    
    figure(f2)
    plot(1:9,mean(rankedRFs_CWnorm_NH{K},1),'k.-')
    
    figure(f3)
    plot(1:9,mean(rankedRFs_sponNorm_NH{K},1),'k.-')
end

for K=1:length(paths_EN)
    figure(f1)
    plot(1:9,mean(rankedRFs_EN{K},1),'r.-')
    
    figure(f2)
    plot(1:9,mean(rankedRFs_CWnorm_EN{K},1),'r.-')
    
    figure(f3)
    plot(1:9,mean(rankedRFs_sponNorm_EN{K},1),'r.-')
end
figure(f1)
title('ranked RFs, raw median')

figure(f2)
title('ranked RFs, norm to CW')

figure(f3)
title('ranked RFs, norm to spont')


rankedRFs_NH=cat(1,rankedRFs_NH{:});
meanRaw_NH=mean(rankedRFs_NH,1);
semRaw_NH=std(rankedRFs_NH,1)/sqrt(size(rankedRFs_NH,1));

rankedRFs_EN=cat(1,rankedRFs_EN{:});
meanRaw_EN=mean(rankedRFs_EN,1);
semRaw_EN=std(rankedRFs_EN,1)/sqrt(size(rankedRFs_EN,1));

figure; hold on
boundedline(1:9,meanRaw_NH,semRaw_NH,'ko-',1:9,meanRaw_EN,semRaw_EN,'ro-','alpha')
title('ranked RFs, raw median')

rankedRFs_CWnorm_NH=cat(1,rankedRFs_CWnorm_NH{:});
meanCWnorm_NH=mean(rankedRFs_CWnorm_NH,1);
semCWnorm_NH=std(rankedRFs_CWnorm_NH,1)/sqrt(size(rankedRFs_CWnorm_NH,1));

rankedRFs_CWnorm_EN=cat(1,rankedRFs_CWnorm_EN{:});
meanCWnorm_EN=mean(rankedRFs_CWnorm_EN,1);
semCWnorm_EN=std(rankedRFs_CWnorm_EN,1)/sqrt(size(rankedRFs_CWnorm_EN,1));

figure; hold on
boundedline(1:9,meanCWnorm_NH,semCWnorm_NH,'ko-',1:9,meanCWnorm_EN,semCWnorm_EN,'ro-','alpha')
title('ranked RFs, norm to CW')

rankedRFs_sponNorm_NH=cat(1,rankedRFs_sponNorm_NH{:});
meansponNorm_NH=mean(rankedRFs_sponNorm_NH,1);
semsponNorm_NH=std(rankedRFs_sponNorm_NH,1)/sqrt(size(rankedRFs_sponNorm_NH,1));

rankedRFs_sponNorm_EN=cat(1,rankedRFs_sponNorm_EN{:});
meansponNorm_EN=mean(rankedRFs_sponNorm_EN,1);
semsponNorm_EN=std(rankedRFs_sponNorm_EN,1)/sqrt(size(rankedRFs_sponNorm_EN,1));

figure; hold on
boundedline(1:9,meansponNorm_NH,semsponNorm_NH,'ko-',1:9,meansponNorm_EN,semsponNorm_EN,'ro-','alpha')
title('ranked RFs, norm to spont')

inds_whisk_NH=repmat(1:9,size(rankedRFs_NH,1),1);
inds_whisk_EN=repmat(1:9,size(rankedRFs_EN,1),1);
EN=cell(size(rankedRFs_EN));
EN(:)={'EN'};
NH=cell(size(rankedRFs_NH));
NH(:)={'NH'};

testV=[rankedRFs_NH(:); rankedRFs_EN(:)];
whiskG=[inds_whisk_NH(:); inds_whisk_EN(:)];
mice=[NH(:); EN(:)];

[P,T,STATS,TERMS]=anovan(testV,{whiskG,mice},'varnames',{'whisker rank','group'},'model','interaction')
pvals.rankedRFraw=P(2);

testV=[rankedRFs_sponNorm_NH(:); rankedRFs_sponNorm_EN(:)];
[P,T,STATS,TERMS]=anovan(testV,{whiskG,mice},'varnames',{'whisker rank','group'},'model','interaction')
pvals.rankedRFnorm=P(2);
%%
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
                        deltaF,sampRate,whiskPref ] = load_NPsub_data_L4( pathNames{J},2 );
                    %                     traceByStim=traceByStim(2);
                    %                     sponTrace=sponTrace(2);
                    %                     permTestResults=permTestResults(2);
                    %                     whiskPref=whiskPref(2);
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


    function [mean_dist,mean_percMatch,sem_percMatch,percMatch_byBarr,percMatch,num_ROIs,CWtunedInds_binned,binEdges]=plot_CWtunedByRad_num(pathNames,npSub,binEdges)
        
        
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
                        deltaF,sampRate,whiskPref ] = load_NPsub_data_L4( pathNames{J},2 );
                    %                     traceByStim=traceByStim(2);
                    %                     sponTrace=sponTrace(2);
                    %                     permTestResults=permTestResults(2);
                    %                     whiskPref=whiskPref(2);
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
            CWtunedInds_binned{i}=CWtunedInds(binInds==i);
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
                        deltaF,sampRate,whiskPref ] = load_NPsub_data_L4( pathNames{J},2 );
                    %                     traceByStim=traceByStim(2);
                    %                     sponTrace=sponTrace(2);
                    %                     permTestResults=permTestResults(2);
                    %                     whiskPref=whiskPref(2);
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
                        deltaF,sampRate,whiskPref ] = load_NPsub_data_L4( pathNames{J},2 );
                    %                     traceByStim=traceByStim(2);
                    %                     sponTrace=sponTrace(2);
                    %                     permTestResults=permTestResults(2);
                    %                     whiskPref=whiskPref(2);
            end
            [percentBWmatch{J},percentBWmatch_num{J}] = find_BWmatchCW( ROIsInBarrel,whiskPref,traceByStim,framesEvoked );
        end
        percentBWmatchCW=horzcat(percentBWmatch{:});
        percentBWmatchCW_num=horzcat(percentBWmatch_num{:});
        
    end

    function [rankedRFs,rankedRFs_CWnorm,rankedRFs_sponNorm]=plot_rankedRFs(pathNames,npSub)
        for J=1:length(pathNames)
            switch npSub
                case 0
                    [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                        dists,ROIsInBarrel,ROItoBarrel,ROI_positions,Stimuli,...
                        deltaF,sampRate,whiskPref ] = load_nonNPsub_data( pathNames{J} );
                case 1
                    [cells,traceByStim,sponTrace,framesEvoked,permTestResults,...
                        dists,ROIsInBarrel,ROItoBarrel,ROI_positions,Stimuli,...
                        deltaF,sampRate,whiskPref ] = load_NPsub_data_L4( pathNames{J},2 );
            end
            ROIs=find_sigROIs(permTestResults,traceByStim);
            [ sponTrace ] = make_sponTrace( Stimuli,sampRate(1),deltaF,0.5,4 );
            [rankedRFs{J},rankedRFs_CWnorm{J},rankedRFs_sponNorm{J}] = find_rankedRFs( traceByStim,sponTrace,ROIs,framesEvoked );
        end
    end
end