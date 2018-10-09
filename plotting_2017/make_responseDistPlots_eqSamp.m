function [results_EN,results_NH]=make_responseDistPlots_eqSamp( paths_EN,paths_NH,npSub,cellsToPlot_EN,cellsToPlot_NH )

% only use fields that have at least 1 ROIs
if ~isempty(cellsToPlot_EN)
    useEN=cellfun(@(x)length(x)>1,cellsToPlot_EN);
    paths_EN=paths_EN(useEN);
    cellsToPlot_EN=cellsToPlot_EN(useEN);
end

if ~isempty(cellsToPlot_NH)
    useNH=cellfun(@(x)length(x)>1,cellsToPlot_NH);
    paths_NH=paths_NH(useNH);
    cellsToPlot_NH=cellsToPlot_NH(useNH);
end

%%

[ blank_rawDF_EN,BW_rawDF_EN,nonBW_rawDF_EN,CW_rawDF_EN,BW_Z_EN,nonBW_Z_EN,...
    CW_Z_EN,numBW_EN,numSigW_EN,fano_by_whisk_EN,medianRs_EN,fano_CW_EN,SW_Z_EN,SW_rawDF_EN]=find_response_dists(paths_EN,npSub,cellsToPlot_EN);
[ blank_rawDF_NH,BW_rawDF_NH,nonBW_rawDF_NH,CW_rawDF_NH,BW_Z_NH,nonBW_Z_NH,...
    CW_Z_NH,numBW_NH,numSigW_NH,fano_by_whisk_NH,medianRs_NH,fano_CW_NH,SW_Z_NH,SW_rawDF_NH]=find_response_dists(paths_NH,npSub,cellsToPlot_NH);

%% return means, sem

results_EN.meanBlank_EN=mean(blank_rawDF_EN);
results_EN.semBlank_EN=std(blank_rawDF_EN)/sqrt(length(blank_rawDF_EN));
results_EN.meanCW_Z_EN=mean(CW_Z_EN);
results_EN.semCW_Z_EN=std(CW_Z_EN)/sqrt(length(SW_Z_EN));
results_EN.meanSW_Z_EN=mean(SW_Z_EN);
results_EN.semSW_Z_EN=std(SW_Z_EN)/sqrt(length(SW_Z_EN));
results_EN.meanBW_Z_EN=mean(BW_Z_EN);
results_EN.semBW_Z_EN=std(BW_Z_EN)/sqrt(length(BW_Z_EN));
results_EN.meanCW_rawDF_EN=mean(CW_rawDF_EN);
results_EN.semCW_rawDF_EN=std(CW_rawDF_EN)/sqrt(length(SW_rawDF_EN));
results_EN.meanSW_rawDF_EN=mean(SW_rawDF_EN);
results_EN.semSW_rawDF_EN=std(SW_rawDF_EN)/sqrt(length(SW_rawDF_EN));
results_EN.meanBW_rawDF_EN=mean(BW_rawDF_EN);
results_EN.semBW_rawDF_EN=std(BW_rawDF_EN)/sqrt(length(BW_rawDF_EN));

results_NH.meanBlank_NH=mean(blank_rawDF_NH);
results_NH.semBlank_NH=std(blank_rawDF_NH)/sqrt(length(blank_rawDF_NH));
results_NH.meanCW_Z_NH=mean(CW_Z_NH);
results_NH.semCW_Z_NH=std(CW_Z_NH)/sqrt(length(SW_Z_NH));
results_NH.meanSW_Z_NH=mean(SW_Z_NH);
results_NH.semSW_Z_NH=std(SW_Z_NH)/sqrt(length(SW_Z_NH));
results_NH.meanBW_Z_NH=mean(BW_Z_NH);
results_NH.semBW_Z_NH=std(BW_Z_NH)/sqrt(length(BW_Z_NH));
results_NH.meanCW_rawDF_NH=mean(CW_rawDF_NH);
results_NH.semCW_rawDF_NH=std(CW_rawDF_NH)/sqrt(length(SW_rawDF_NH));
results_NH.meanSW_rawDF_NH=mean(SW_rawDF_NH);
results_NH.semSW_rawDF_NH=std(SW_rawDF_NH)/sqrt(length(SW_rawDF_NH));
results_NH.meanBW_rawDF_NH=mean(BW_rawDF_NH);
results_NH.semBW_rawDF_NH=std(BW_rawDF_NH)/sqrt(length(BW_rawDF_NH));

%% box plots

blankData{1}=blank_rawDF_NH;
blankData{2}=blank_rawDF_EN;

figure; hold on
plotSpread(blankData,'showMM',4,'distributionColors',{'k','r'},'xNames',{'Ct, blank','En, blank'})
ylabel('median dF/F')

CW_Z{1}=CW_Z_NH;
CW_Z{2}=CW_Z_EN;

SW_Z{1}=SW_Z_NH;
SW_Z{2}=SW_Z_EN;

BW_Z{1}=BW_Z_NH;
BW_Z{2}=BW_Z_EN;

figure; hold on
CW=plotSpread(CW_Z,'showMM',4,'distributionColors',{'k','r'},'xNames',{'Ct, CW','En, CW'},'xValues',[1 1.25],'spreadWidth',0.25);
SW=plotSpread(SW_Z,'showMM',4,'distributionColors',{'k','r'},'xNames',{'Ct, SW','En, SW'},'xValues',[2 2.25],'spreadWidth',0.25)
BW=plotSpread(BW_Z,'showMM',4,'distributionColors',{'k','r'},'xNames',{'Ct, BW','En, BW'},'xValues',[3 3.25],'spreadWidth',0.25)
CW{2}(1).Color='g';
CW{2}(2).Color='g';
CW{2}(1).LineWidth=2.5;
CW{2}(2).LineWidth=2.5;
SW{2}(1).Color='g';
SW{2}(2).Color='g';
SW{2}(1).LineWidth=2.5;
SW{2}(2).LineWidth=2.5;
BW{2}(1).Color='g';
BW{2}(2).Color='g';
BW{2}(1).LineWidth=2.5;
BW{2}(2).LineWidth=2.5;

ax=gca;
ax.XTick=[1 1.25 2 2.25 3 3.25];
ax.XTickLabel={'Ct, CW','En, CW','Ct, SW','En, SW','Ct, BW','En, BW'};
ylabel('Z-scored median dF/F')

CW_rawDF{1}=CW_rawDF_NH;
CW_rawDF{2}=CW_rawDF_EN;

SW_rawDF{1}=SW_rawDF_NH;
SW_rawDF{2}=SW_rawDF_EN;

BW_rawDF{1}=BW_rawDF_NH;
BW_rawDF{2}=BW_rawDF_EN;

figure; hold on
CW=plotSpread(CW_rawDF,'showMM',4,'distributionColors',{'k','r'},'xNames',{'Ct, CW','En, CW'},'xValues',[1 1.25],'spreadWidth',0.25)
SW=plotSpread(SW_rawDF,'showMM',4,'distributionColors',{'k','r'},'xNames',{'Ct, SW','En, SW'},'xValues',[2 2.25],'spreadWidth',0.25)
BW=plotSpread(BW_rawDF,'showMM',4,'distributionColors',{'k','r'},'xNames',{'Ct, BW','En, BW'},'xValues',[3 3.25],'spreadWidth',0.25)
CW{2}(1).Color='g';
CW{2}(2).Color='g';
CW{2}(1).LineWidth=2.5;
CW{2}(2).LineWidth=2.5;
SW{2}(1).Color='g';
SW{2}(2).Color='g';
SW{2}(1).LineWidth=2.5;
SW{2}(2).LineWidth=2.5;
BW{2}(1).Color='g';
BW{2}(2).Color='g';
BW{2}(1).LineWidth=2.5;
BW{2}(2).LineWidth=2.5;

ax=gca;
ax.XTick=[1 1.25 2 2.25 3 3.25];
ax.XTickLabel={'Ct, CW','En, CW','Ct, SW','En, SW','Ct, BW','En, BW'};
ylabel('median dF/F')

%% fano factor plots

fano_CW_EN=cat(2,fano_CW_EN{:});
fano_CW_NH=cat(2,fano_CW_NH{:});

figure; hold on
plot_2cdfs(fano_CW_NH,fano_CW_EN);
xlabel('fano factor, CW')

fano_all_EN=cat(1,fano_by_whisk_EN{:});
fano_all_NH=cat(1,fano_by_whisk_NH{:});

fano_blank_EN=fano_all_EN(:,10);
fano_blank_NH=fano_all_NH(:,10);

fano_all_EN=fano_all_EN(:,1:9);
fano_all_NH=fano_all_NH(:,1:9);

figure; hold on
plot_2cdfs(fano_blank_NH,fano_blank_EN);
xlabel('fano factor, blank trials');

figure; hold on
plot_2cdfs(fano_all_NH(:),fano_all_EN(:));
xlabel('fano factor, all trials');

fano_field_EN=cellfun(@(x)mean(mean(x(:,1:9))),fano_by_whisk_EN);
fano_field_NH=cellfun(@(x)mean(mean(x(:,1:9))),fano_by_whisk_NH);
notBoxPlot_ENvNH( fano_field_EN,fano_field_NH )
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

%% Z-scored dists

figure; hold on
subplot(1,3,1); hold on
p1=cdfplot(CW_Z_NH);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(CW_Z_EN)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
ylabel('')
xlabel('Z-Scored median dF/F')
pval=permutationTest(CW_Z_EN,CW_Z_NH,10000);
title(strcat('CW response,p=',num2str(pval)))

subplot(1,3,2); hold on
p1=cdfplot(SW_Z_NH);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(SW_Z_EN)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
ylabel('')
xlabel('Z-Scored median dF/F')
pval=permutationTest(SW_Z_EN,SW_Z_NH,10000);
title(strcat('mean SW response,p=',num2str(pval)))



subplot(1,3,3); hold on
p2=cdfplot(BW_Z_NH);
p1=cdfplot(BW_Z_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(BW_Z_EN,BW_Z_NH,100000);
title(strcat('BW response, p=',num2str(pval)))

%%
figure; hold on
subplot(1,3,1); hold on
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


subplot(1,3,2); hold on
p1=cdfplot(SW_rawDF_NH);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(SW_rawDF_EN)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
legend('NH','EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(SW_rawDF_EN,SW_rawDF_NH,10000);
title(strcat('mean SW response,p=',num2str(pval)))

subplot(1,3,3); hold on
p2=cdfplot(BW_rawDF_NH);
p1=cdfplot(BW_rawDF_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(BW_rawDF_EN,BW_rawDF_NH,10000);
title(strcat('BW response, p=',num2str(pval)))

%%
% tmp=gca;
% tmp.XLim(2)=0.15
% 
figure; hold on
p2=cdfplot(nonBW_rawDF_NH);
p1=cdfplot(nonBW_rawDF_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(nonBW_rawDF_EN,nonBW_rawDF_NH,10000);
title(strcat('non-BW response, p=',num2str(pval)))


%%


figure; hold on
p2=cdfplot(nonBW_Z_NH);
p1=cdfplot(nonBW_Z_EN);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('NH','EN')
ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(nonBW_Z_EN,nonBW_Z_NH,10000);
title(strcat('non-BW response, p=',num2str(pval)))

%%



%% SW responses


%% median response analysis for Enriched animals


    





end

