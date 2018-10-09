function make_responseDistPlots( paths_EN,paths_NH,npSub,rval,cellsToPlot_EN,cellsToPlot_NH )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[ blank_rawDF_EN,BW_rawDF_EN,bestSWs_rawDF_EN,worstSWs_rawDF_EN,CW_rawDF_EN,...
    BW_Z_EN,bestSWs_Z_EN,worstSWs_Z_EN,CW_Z_EN,numBW_EN,numSigW_EN,fano_by_whisk_EN,...
    medianRs_EN,fano_CW_EN]=find_response_dists(paths_EN,npSub,rval,cellsToPlot_EN);

[ blank_rawDF_NH,BW_rawDF_NH,bestSWs_rawDF_NH,worstSWs_rawDF_NH,CW_rawDF_NH,...
    BW_Z_NH,bestSWs_Z_NH,worstSWs_Z_NH,CW_Z_NH,numBW_NH,numSigW_NH,fano_by_whisk_NH,...
    medianRs_NH,fano_CW_NH]=find_response_dists(paths_NH,npSub,rval,cellsToPlot_NH);


figure; hold on
plot_scatterRLine(worstSWs_rawDF_EN,bestSWs_rawDF_EN);
% alpha_scatter(worstSWs_rawDF_NH,bestSWs_rawDF_NH,'markerSize,',8,'markerColor',[0 0 0]);
% alpha_scatter(worstSWs_rawDF_EN,bestSWs_rawDF_EN,'markerSize,',8,'markerColor',[1 0 0]);
% [r2_EN,p_EN,slope_EN]=plot_RLine_only(worstSWs_rawDF_EN,bestSWs_rawDF_EN,'r');
% [r2_NH,p_NH,slope_NH]=plot_RLine_only(worstSWs_rawDF_NH,bestSWs_rawDF_NH,'k');
xlabel('EN worst SWs, raw')
ylabel('EN best SWs,raw')


[EN,~,binEdgesX,binEdgesY] = scatterPlot_marginals( worstSWs_rawDF_EN,bestSWs_rawDF_EN,[],[] );
NH=scatterPlot_marginals(worstSWs_rawDF_NH,bestSWs_rawDF_NH,binEdgesX,binEdgesY);

ENscatt=EN.Children(3);
NHscatt=NH.Children(3);

linkaxes([ENscatt,NHscatt],'xy')

p_Worst=permutationTest(worstSWs_rawDF_EN,worstSWs_rawDF_NH,10000)
p_Best=permutationTest(bestSWs_rawDF_EN,bestSWs_rawDF_NH,10000)

% test for difference in regression coefficients
ln_EN=fitlm(worstSWs_rawDF_EN,bestSWs_rawDF_EN)
ln_NH=fitlm(worstSWs_rawDF_NH,bestSWs_rawDF_NH)

%compute tstat
T=abs(ln_EN.Coefficients.Estimate(2)-ln_NH.Coefficients.Estimate(2))/sqrt((ln_EN.Coefficients.SE(2)^2)+(ln_NH.Coefficients.SE(2)^2))
DFE=2*(ln_EN.DFE+ln_NH.DFE)-2
%p-value
p = 1-tcdf(T,DFE);

figure; hold on
plot_scatterRLine(worstSWs_rawDF_NH,bestSWs_rawDF_NH);
xlabel('NH worst SWs, raw')
ylabel('NH best SWs, raw')

bestWorstIdx_EN=bestSWs_rawDF_EN./worstSWs_rawDF_EN;
bestWorstIdx_NH=bestSWs_rawDF_NH./worstSWs_rawDF_NH;
figure; hold on
plot_2cdfs(bestWorstIdx_NH,bestWorstIdx_EN);

figure; hold on
plot_scatterRLine(worstSWs_Z_EN,bestSWs_Z_EN);
xlabel('EN worst SWs, Z')
ylabel('EN best SWs, Z')

figure; hold on
plot_scatterRLine(worstSWs_Z_NH,bestSWs_Z_NH);
xlabel('NH worst SWs, Z')
ylabel('NH best SWs, Z')


figure; hold on
plot_scatterRLine(blank_rawDF_EN,bestSWs_rawDF_EN);
xlabel('EN spont')
ylabel('EN best SWs, raw')

figure; hold on
plot_scatterRLine(blank_rawDF_NH,bestSWs_rawDF_NH);
xlabel('NH spont')
ylabel('NH best SWs, raw')

figure; hold on
plot_scatterRLine(blank_rawDF_EN,worstSWs_rawDF_EN);
xlabel('EN spont')
ylabel('EN worst SWs, raw')

figure; hold on
plot_scatterRLine(blank_rawDF_NH,worstSWs_rawDF_NH);
xlabel('NH spont')
ylabel('NH worst SWs, raw')

figure; hold on
plot_scatterRLine(blank_rawDF_EN,bestSWs_Z_EN);
xlabel('EN spont')
ylabel('EN best SWs, Z')

figure; hold on
plot_scatterRLine(blank_rawDF_NH,bestSWs_Z_NH);
xlabel('NH spont')
ylabel('NH best SWs, Z')

figure; hold on
plot_scatterRLine(blank_rawDF_EN,worstSWs_Z_EN);
xlabel('EN spont')
ylabel('EN worst SWs, Z')

figure; hold on
plot_scatterRLine(blank_rawDF_NH,worstSWs_Z_NH);
xlabel('NH spont')
ylabel('NH worst SWs, Z')

%%
figure; hold on
plot_scatterRLine(worstSWs_rawDF_EN,BW_rawDF_EN);
xlabel('EN worst SWs, raw')
ylabel('EN BW,raw')

figure; hold on
plot_scatterRLine(worstSWs_rawDF_NH,BW_rawDF_NH);
xlabel('NH worst SWs, raw')
ylabel('NH BW, raw')

figure; hold on
plot_scatterRLine(bestSWs_rawDF_EN,BW_rawDF_EN);
xlabel('EN best SWs, raw')
ylabel('EN BW,raw')

figure; hold on
plot_scatterRLine(bestSWs_rawDF_NH,BW_rawDF_NH);
xlabel('NH best SWs, raw')
ylabel('NH BW, raw')

% numSigW=[numSigW_EN; numSigW_NH];
% numBW=[numBW_EN; numBW_NH];
% BW_rawDF=[BW_rawDF_EN, BW_rawDF_NH];
% 
% 
% for i=1:9
%     meanBWr_numSigW(i)=mean(BW_rawDF(numSigW==i));
%     semBWr_numSigW(i)=std(BW_rawDF(numSigW==i))/sqrt(sum(numSigW==i));
%     meanBWr_numBW(i)=mean(BW_rawDF(numBW==i));
%     semBWr_numBW(i)=std(BW_rawDF(numBW==i))/sqrt(sum(numBW==i));
%     BWr_numSigW{i}=BW_rawDF(numSigW==i);
%      BWr_numBW{i}=BW_rawDF(numBW==i);
%      
%      % for NH only
%      meanBWr_numSigW_NH(i)=mean(BW_rawDF_NH(numSigW_NH==i));
%     semBWr_numSigW_NH(i)=std(BW_rawDF_NH(numSigW_NH==i))/sqrt(sum(numSigW_NH==i));
%     meanBWr_numBW_NH(i)=mean(BW_rawDF_NH(numBW_NH==i));
%     semBWr_numBW_NH(i)=std(BW_rawDF_NH(numBW_NH==i))/sqrt(sum(numBW_NH==i));
%     BWr_numSigW_NH{i}=BW_rawDF_NH(numSigW_NH==i);
%      BWr_numBW_NH{i}=BW_rawDF_NH(numBW_NH==i);
% 
%      % for EN only
%      meanBWr_numSigW_EN(i)=mean(BW_rawDF_EN(numSigW_EN==i));
%     semBWr_numSigW_EN(i)=std(BW_rawDF_EN(numSigW_EN==i))/sqrt(sum(numSigW_EN==i));
%     meanBWr_numBW_EN(i)=mean(BW_rawDF_EN(numBW_EN==i));
%     semBWr_numBW_EN(i)=std(BW_rawDF_EN(numBW_EN==i))/sqrt(sum(numBW_EN==i));
%     BWr_numSigW_EN{i}=BW_rawDF_EN(numSigW_EN==i);
%      BWr_numBW_EN{i}=BW_rawDF_EN(numBW_EN==i);
% 
% end
% % 
% % figure; hold on
% % plotSpread(BWr_numSigW,'distributionColors','k');
% % boundedline(1:9,meanBWr_numSigW,semBWr_numSigW,'ro-','alpha')
% % xlabel('# whiskers in RF')
% % ylabel('median response to BW')
% % [p_sigW,tbl,stats_sigW] = anova1(BW_rawDF,numSigW);
% % title(['p= ',num2str(p_sigW)])
% % 
% % figure; hold on
% % plotSpread(BWr_numBW,'distributionColors','k');
% % boundedline(1:9,meanBWr_numBW,semBWr_numBW,'ro-','alpha')
% % xlabel('# eq BWs')
% % ylabel('median response to BW')
% % [p_numBW,tbl,stats_numBW] = anova1(BW_rawDF,numBW);
% % title(['p= ',num2str(p_numBW)])
% % tmp=multcompare(stats_numBW)
% 
% 
% figure; hold on
% plotSpread(BWr_numSigW_NH,'distributionColors','k');
% plotSpread(BWr_numSigW_EN,'distributionColors','r');
% boundedline(1:9,meanBWr_numSigW_NH,semBWr_numSigW_NH,'ko-','alpha')
% boundedline(1:9,meanBWr_numSigW_EN,semBWr_numSigW_EN,'ro-','alpha')
% xlabel('# whiskers in RF')
% ylabel('median response to BW')
% [p_sigW_NH,tbl,stats_numSigW_NH] = anova1(BW_rawDF_NH,numSigW_NH);
% [p_sigW_EN,tbl,stats_numSigW_EN] = anova1(BW_rawDF_EN,numSigW_EN);
% 
% 
% % [p_sigW,tbl,stats_sigW] = anova1(BW_rawDF,numSigW);
% % title(['p= ',num2str(p_sigW)])
% 
% figure; hold on
% plotSpread(BWr_numBW_NH,'distributionColors','k');
% plotSpread(BWr_numBW_EN,'distributionColors','r');
% boundedline(1:8,meanBWr_numBW_NH(1:8),semBWr_numBW_NH(1:8),'ko-','alpha')
% boundedline(1:9,meanBWr_numBW_EN,semBWr_numBW_EN,'ro-','alpha')
% xlabel('# eq BWs')
% ylabel('median response to BW')
% [p_numBW,tbl,stats_numBW_NH] = anova1(BW_rawDF_NH,numBW_NH);
% [p_numBW,tbl,stats_numBW_EN] = anova1(BW_rawDF_EN,numBW_EN);
% 
% 
% % title(['p= ',num2str(p_numBW)])
% multcompare(stats_numSigW_NH)


% figure; hold on
% plot(numSigW,BW_rawDF,'ko')
% plot(1:9,meanBWr_numSigW,'ro-')
% xlabel('# whiskers driving sig. response')
% ylabel('median response to BW')
% 
% figure; hold on
% plot(numBW,BW_rawDF,'ko')
% plot(1:9,meanBWr_numBW,'ro-')
% xlabel('# equivalent BWs')
% ylabel('median response to BW')
% %% fano factor plots
% 
% fano_CW_EN=cat(2,fano_CW_EN{:});
% fano_CW_NH=cat(2,fano_CW_NH{:});
% 
% figure; hold on
% plot_2cdfs(fano_CW_NH,fano_CW_EN);
% xlabel('fano factor, CW')
% 
% fano_all_EN=cat(1,fano_by_whisk_EN{:});
% fano_all_NH=cat(1,fano_by_whisk_NH{:});
% 
% fano_blank_EN=fano_all_EN(:,10);
% fano_blank_NH=fano_all_NH(:,10);
% 
% fano_all_EN=fano_all_EN(:,1:9);
% fano_all_NH=fano_all_NH(:,1:9);
% 
% figure; hold on
% plot_2cdfs(fano_blank_NH,fano_blank_EN);
% xlabel('fano factor, blank trials');
% 
% figure; hold on
% plot_2cdfs(fano_all_NH(:),fano_all_EN(:));
% xlabel('fano factor, all trials');
% 
% fano_field_EN=cellfun(@(x)mean(mean(x(:,1:9))),fano_by_whisk_EN);
% fano_field_NH=cellfun(@(x)mean(mean(x(:,1:9))),fano_by_whisk_NH);
% notBoxPlot_ENvNH( fano_field_EN,fano_field_NH )
%%
figure; hold on
p2=cdfplot(blank_rawDF_NH);
p1=cdfplot(blank_rawDF_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(blank_rawDF_EN,blank_rawDF_NH,10000);
title(strcat('spontaneous dF/F, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.2



figure; hold on
p2=cdfplot(BW_rawDF_NH);
p1=cdfplot(BW_rawDF_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
% tmp=gca;
% tmp.XLim(2)=0.2;



ylabel('')
xlabel('median dF/F')
pval=permutationTest(BW_rawDF_EN,BW_rawDF_NH,10000);
title(strcat('BW response, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.15
% 
figure; hold on
p2=cdfplot(bestSWs_rawDF_NH);
p1=cdfplot(bestSWs_rawDF_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(bestSWs_rawDF_EN,bestSWs_rawDF_NH,10000);
title(strcat('responses to best SWs, p=',num2str(pval)))

figure; hold on
p2=cdfplot(worstSWs_rawDF_NH);
p1=cdfplot(worstSWs_rawDF_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(worstSWs_rawDF_EN,worstSWs_rawDF_NH,10000);
title(strcat('responses to worst SWs, p=',num2str(pval)))



%%
figure; hold on
p2=cdfplot(BW_Z_NH);
p1=cdfplot(BW_Z_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
% tmp=gca;
% tmp.XLim(2)=0.2;



ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(BW_Z_EN,BW_Z_NH,100000);
title(strcat('BW response, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.15

figure; hold on
p2=cdfplot(bestSWs_Z_NH);
p1=cdfplot(bestSWs_Z_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(bestSWs_Z_EN,bestSWs_Z_NH,10000);
title(strcat('responses to best SWs, p=',num2str(pval)))

figure; hold on
p2=cdfplot(worstSWs_Z_NH);
p1=cdfplot(worstSWs_Z_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(worstSWs_Z_EN,worstSWs_Z_NH,10000);
title(strcat('responses to worst SWs, p=',num2str(pval)))

%%

figure; hold on
p1=cdfplot(CW_Z_NH);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(CW_Z_EN)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
legend('NH','EN')
ylabel('')
xlabel('Z-Scored median dF/F')
pval=permutationTest(CW_Z_EN,CW_Z_NH,10000);
title(strcat('CW response,p=',num2str(pval)))

figure; hold on
p1=cdfplot(CW_rawDF_NH);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(CW_rawDF_EN)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
legend('NH','EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(CW_rawDF_EN,CW_rawDF_NH,10000);
title(strcat('CW response,p=',num2str(pval)))
%% median response analysis for Enriched animals


%     function [ blank_rawDF,BW_rawDF,nonBW_rawDF,CW_rawDF,BW_Z,nonBW_Z,CW_Z,numBW,numSigW,fano_by_whisk,medianRs,fano_CW]=find_response_dists(pathNames,npSub,rval,cellsToPlot)
%         
%        
%         
%         
%         for K=1:length(pathNames)
%             switch npSub
%                 case 1
%                     [traceByStim,sponTrace,framesEvoked,permTestResults,~,ROIsInBarrel,~,~,Stimuli,deltaF,sampRate,whiskPref ] = load_NPsub_data_L23( pathNames{K},rval );
%                
%                 case 0
%                     [~,traceByStim,sponTrace,framesEvoked,permTestResults,~,ROIsInBarrel,~,~,Stimuli,deltaF,sampRate,whiskPref ] = load_nonNPsub_data( pathNames{K} );
%             end
%             [ sponTraceRaw ] = make_sponTrace_noBLS( Stimuli,sampRate(1),deltaF,0.5,3 );
%             
%             if ~isempty(cellsToPlot)
%                 sigCells=cellsToPlot{K};
%             else
%                 sigCells=find_sigROIs( permTestResults,traceByStim ); %only include significantly responsive cells
%             end
%             cellNames=fieldnames(traceByStim);
%             [ rawDF(K),Zscored(K) ] = activity_medians(sigCells,traceByStim,sponTraceRaw,sponTrace,framesEvoked );
%             
%             [ ~,Zscored_CW{K},medians_raw_CW{K},fano_CW{K} ] = find_CW_responses_Z(sigCells,ROIsInBarrel,traceByStim,sponTrace,framesEvoked);
%             numBW{K}=cellfun(@(x)length(whiskPref.(x){1}),sigCells,'Uni',1);
%             numSigW{K}=cellfun(@(x)length(horzcat(whiskPref.(x){:})),sigCells,'Uni',1);
%             
%            % fano factor
%            [ fano_by_whisk{K},medianRs{K} ] = fano_factor( sigCells,traceByStim,sponTrace,framesEvoked );
%            
% %            figure; hold on
% %            plot_scatterRLine(medianRs{K}(:),fano_by_whisk{K}(:));
% %            xlabel('median deltaF/F')
% %            ylabel('fano factor')
% %            
%             clear traceByStim
%             clear sponTrace
%             clear sigCells
%             
%             
%         end
%         blank_rawDF=cat(2,rawDF(:).blank);
%         BW_rawDF=cat(2,rawDF(:).BWr);
%         nonBW_rawDF=cat(2,rawDF(:).nonBW);
%         CW_rawDF=cat(2,medians_raw_CW{:});
%         
%         BW_Z=cat(2,Zscored(:).BWr);
%         nonBW_Z=cat(2,Zscored(:).nonBW);
%         CW_Z=cat(2,Zscored_CW{:});
%         
%         numBW=cat(1,numBW{:});
%         numSigW=cat(1,numSigW{:});
%          
%         % fano factor
%     end






end

