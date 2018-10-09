function pvals=plot_ENvsNH_OWvsMW(paths_EN,paths_NH,npSub )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% [OWtuned_res,MWtuned_res]=plot_MWtunedDists( pathNames,npSub )
[CWtuned_EN,SWtuned_EN,MWtuned_CW_EN,OWtuned_CW_EN,propMW_EN]=plot_MWtunedDists( paths_EN,npSub );
[CWtuned_NH,SWtuned_NH,MWtuned_CW_NH,OWtuned_CW_NH,propMW_NH]=plot_MWtunedDists( paths_NH,npSub );

%%
[ pval_MW] = notBoxPlot_ENvNH( MWtuned_CW_EN,MWtuned_CW_NH )
title(['proportion of multi-whisker cells tuned to CW, p=',num2str(pval_MW)])
[ pval_OW] = notBoxPlot_ENvNH( OWtuned_CW_EN,OWtuned_CW_NH )
title(['proportion of single-whisker cells tuned to CW, p=',num2str(pval_OW)])
[ pval_propMW] = notBoxPlot_ENvNH( propMW_EN,propMW_NH )
title(['proportion of multi-whisker cells in field, p=',num2str(pval_propMW)])

figure;
hold on
plot_4cdfs(SWtuned_NH.D2dist,SWtuned_EN.D2dist,CWtuned_NH.D2dist,CWtuned_EN.D2dist)
legend('MW tuned NH','MW tuned EN','OW tuned NH','OW tuned EN')
ylabel('')
xlabel('um')
title('distance to D2')

%% compare CW-tuned ROIs between EN and NH
figure;
hold on
p2=cdfplot(CWtuned_NH.blank_rawDF);
p1=cdfplot(CWtuned_EN.blank_rawDF);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('CW tuned NH','CW tuned EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(CWtuned_NH.blank_rawDF,CWtuned_EN.blank_rawDF,10000);
title(strcat('spontaneous dF/F, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.2



figure;
subplot(2,3,1);
hold on
p2=cdfplot(CWtuned_NH.BW_rawDF);
p1=cdfplot(CWtuned_EN.BW_rawDF);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';

% tmp=gca;
% tmp.XLim(2)=0.2;



ylabel('')
xlabel('median dF/F')
pval=permutationTest(CWtuned_NH.BW_rawDF,CWtuned_EN.BW_rawDF,100000);
title(strcat('BW response, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.15
% 

subplot(2,3,2);
hold on

p2=cdfplot(CWtuned_NH.nonBW_rawDF);
p1=cdfplot(CWtuned_EN.nonBW_rawDF);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
ylabel('')
xlabel('median dF/F')
pval=permutationTest(CWtuned_NH.nonBW_rawDF,CWtuned_EN.nonBW_rawDF,10000);
title(strcat('non-BW response, p=',num2str(pval)))



subplot(2,3,4);
hold on
p2=cdfplot(CWtuned_NH.BW_Z);
p1=cdfplot(CWtuned_EN.BW_Z);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';

ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(CWtuned_NH.BW_Z,CWtuned_EN.BW_Z,100000);
title(strcat('BW response, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.15

subplot(2,3,5);
hold on
p2=cdfplot(CWtuned_NH.nonBW_Z);
p1=cdfplot(CWtuned_EN.nonBW_Z);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('CW tuned NH','CW tuned EN')
ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(CWtuned_NH.nonBW_Z,CWtuned_EN.nonBW_Z,10000);
title(strcat('non-BW response, p=',num2str(pval)))



subplot(2,3,6);
hold on
p1=cdfplot(CWtuned_NH.CW_Z);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(CWtuned_EN.CW_Z)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
legend('CW tuned NH','CW tuned EN')
ylabel('')
xlabel('Z-Scored median dF/F')
pval=permutationTest(CWtuned_NH.CW_Z,CWtuned_EN.CW_Z,10000);
title(strcat('CW response,p=',num2str(pval)))

subplot(2,3,3);
hold on
p1=cdfplot(CWtuned_NH.CW_rawDF);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(CWtuned_EN.CW_rawDF)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
ylabel('')
xlabel('median dF/F')
pval=permutationTest(CWtuned_NH.CW_rawDF,CWtuned_EN.CW_rawDF,10000);
title(strcat('CW response,p=',num2str(pval)))


%% compare SW-tuned ROIs between EN and NH
figure;
hold on
p2=cdfplot(SWtuned_NH.blank_rawDF);
p1=cdfplot(SWtuned_EN.blank_rawDF);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('SW tuned NH','SW tuned EN')
ylabel('')
xlabel('median dF/F')
pval=permutationTest(SWtuned_NH.blank_rawDF,SWtuned_EN.blank_rawDF,10000);
title(strcat('spontaneous dF/F, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.2



figure;
subplot(2,3,1);
hold on
p2=cdfplot(SWtuned_NH.BW_rawDF);
p1=cdfplot(SWtuned_EN.BW_rawDF);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';

% tmp=gca;
% tmp.XLim(2)=0.2;



ylabel('')
xlabel('median dF/F')
pval=permutationTest(SWtuned_NH.BW_rawDF,SWtuned_EN.BW_rawDF,100000);
title(strcat('BW response, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.15
% 

subplot(2,3,2);
hold on

p2=cdfplot(SWtuned_NH.nonBW_rawDF);
p1=cdfplot(SWtuned_EN.nonBW_rawDF);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
ylabel('')
xlabel('median dF/F')
pval=permutationTest(SWtuned_NH.nonBW_rawDF,SWtuned_EN.nonBW_rawDF,10000);
title(strcat('non-BW response, p=',num2str(pval)))



subplot(2,3,4);
hold on
p2=cdfplot(SWtuned_NH.BW_Z);
p1=cdfplot(SWtuned_EN.BW_Z);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';

ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(SWtuned_NH.BW_Z,SWtuned_EN.BW_Z,100000);
title(strcat('BW response, p=',num2str(pval)))
% tmp=gca;
% tmp.XLim(2)=0.15

subplot(2,3,5);
hold on
p2=cdfplot(SWtuned_NH.nonBW_Z);
p1=cdfplot(SWtuned_EN.nonBW_Z);
p2.LineWidth=2.5;
p2.Color='k';
p1.LineWidth=2.5;
p1.Color='r';
legend('SW tuned NH','SW tuned EN')
ylabel('')
xlabel('Z-scored median dF/F')
pval=permutationTest(SWtuned_NH.nonBW_Z,SWtuned_EN.nonBW_Z,10000);
title(strcat('non-BW response, p=',num2str(pval)))



subplot(2,3,6);
hold on
p1=cdfplot(SWtuned_NH.CW_Z);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(SWtuned_EN.CW_Z)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
legend('CW tuned NH','CW tuned EN')
ylabel('')
xlabel('Z-Scored median dF/F')
pval=permutationTest(SWtuned_NH.CW_Z,SWtuned_EN.CW_Z,10000);
title(strcat('CW response,p=',num2str(pval)))

subplot(2,3,3);
hold on
p1=cdfplot(SWtuned_NH.CW_rawDF);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(SWtuned_EN.CW_rawDF)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
ylabel('')
xlabel('median dF/F')
pval=permutationTest(SWtuned_NH.CW_rawDF,SWtuned_EN.CW_rawDF,10000);
title(strcat('CW response,p=',num2str(pval)))

%% compare CW responses between EN and NH, separately for SW and CW -tuned ROIs
figure; 
subplot(1,2,1);
hold on
p1=cdfplot(CWtuned_NH.CW_Z);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(CWtuned_EN.CW_Z)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
legend('CW tuned NH','CW tuned EN')
ylabel('')
xlabel('Z-Scored median dF/F')
pval=permutationTest(CWtuned_NH.CW_Z,CWtuned_EN.CW_Z,10000);
title(strcat('CW response,p=',num2str(pval)))

subplot(1,2,2)
hold on
p1=cdfplot(SWtuned_NH.CW_Z);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(SWtuned_EN.CW_Z)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
legend('SW tuned NH','SW tuned EN')
ylabel('')
xlabel('Z-Scored median dF/F')
pval=permutationTest(SWtuned_NH.CW_Z,SWtuned_EN.CW_Z,10000);
title(strcat('CW response,p=',num2str(pval)))

figure;
subplot(1,2,1);
hold on
p1=cdfplot(CWtuned_NH.CW_rawDF);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(CWtuned_EN.CW_rawDF)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
ylabel('')
xlabel('median dF/F')
pval=permutationTest(CWtuned_NH.CW_rawDF,CWtuned_EN.CW_rawDF,10000);
title(strcat('CW response,p=',num2str(pval)))
legend('CW tuned NH','CW tuned EN')


subplot(1,2,2);
hold on
p1=cdfplot(SWtuned_NH.CW_rawDF);
p1.LineWidth=2.5;
p1.Color='k'
p2=cdfplot(SWtuned_EN.CW_rawDF)%,'FaceColor','none','EdgeColor','c','LineWidth',2)
p2.LineWidth=2.5
p2.Color='r'
ylabel('')
xlabel('median dF/F')
pval=permutationTest(SWtuned_NH.CW_rawDF,SWtuned_EN.CW_rawDF,10000);
title(strcat('CW response,p=',num2str(pval)))
legend('SW tuned NH','SW tuned EN')

%% RF size plots
figure; 
subplot(2,2,1)
hold on
histogram(CWtuned_NH.numBW,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.4,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(CWtuned_EN.numBW,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
ax=gca;
ax.XTick=1:9;
ax.XTickLabel=1:9;
% ax.Children(1).FaceColor='none';
% ax.Children(2).FaceColor='none';



numBWcounts_NH=histcounts(CWtuned_NH.numBW,1:10);
numBWcounts_EN=histcounts(CWtuned_EN.numBW,1:10);

[h,p,ch]=chisquare_binned(1:9,numBWcounts_NH,numBWcounts_EN)
% [tbl,chi2stat,pval]=crosstab(numBWcounts_NH,numBWcounts_EN);
% [~,p]=ttest2(CWtuned_NH.numBW,CWtuned_EN.numBW);
title(strcat('CW tuned, p=',num2str(p)))
xlabel('# equivalent BWs')
ylabel('fraction of ROIs')


subplot(2,2,3); hold on
histogram(CWtuned_NH.numSigW,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.4,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(CWtuned_EN.numSigW,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
ax=gca;
ax.XTick=1:9;
ax.XTickLabel=1:9;

xlabel('# whiskers in RF')
ylabel('fraction of ROIs')

numSigcounts_NH=histcounts(CWtuned_NH.numSigW,1:10);
numSigcounts_EN=histcounts(CWtuned_EN.numSigW,1:10);

[h,p,ch]=chisquare_binned(1:9,numSigcounts_NH,numSigcounts_EN)
% [tbl,chi2stat,pval]=crosstab(numSigcounts_NH,numSigcounts_EN);

% [~,p]=ttest2(CWtuned_NH.numSigW,CWtuned_EN.numSigW);
title(['CW tuned, p=',num2str(p)])

subplot(2,2,2)
hold on
histogram(SWtuned_NH.numBW,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.4,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(SWtuned_EN.numBW,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');

numBWcounts_NH=histcounts(SWtuned_NH.numBW,1:10);
numBWcounts_EN=histcounts(SWtuned_EN.numBW,1:10);

% [tbl,chi2stat,pval]=crosstab(numBWcounts_NH,numBWcounts_EN);
[h,p,ch]=chisquare_binned(1:9,numBWcounts_NH,numBWcounts_EN)
ax=gca;
ax.XTick=1:9;
ax.XTickLabel=1:9;
% [~,p]=ttest2(SWtuned_NH.numBW,SWtuned_EN.numBW);
title(strcat('SW tuned, p=',num2str(p)))


xlabel('# equivalent BWs')
ylabel('fraction of ROIs')


subplot(2,2,4); hold on
histogram(SWtuned_NH.numSigW,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.4,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(SWtuned_EN.numSigW,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');

numSigcounts_NH=histcounts(SWtuned_NH.numSigW,1:10);
numSigcounts_EN=histcounts(SWtuned_EN.numSigW,1:10);

% [tbl,chi2stat,pval]=crosstab(numSigcounts_NH,numSigcounts_EN);
[h,p,ch]=chisquare_binned(1:9,numSigcounts_NH,numSigcounts_EN)
ax=gca;
ax.XTick=1:9;
ax.XTickLabel=1:9;
xlabel('# whiskers in RF')
ylabel('fraction of ROIs')

% [~,p]=ttest2(SWtuned_NH.numSigW,SWtuned_EN.numSigW);
title(['SW tuned, p=',num2str(p)])
legend('NH','EN')

%% plot all groups together

% spontaneous activity
figure;
hold on
plot_4cdfs(SWtuned_NH.blank_rawDF,SWtuned_EN.blank_rawDF,CWtuned_NH.blank_rawDF,CWtuned_EN.blank_rawDF)
legend('SW tuned NH','SW tuned EN','CW tuned NH','CW tuned EN')
ylabel('')
xlabel('median dF/F')
title('spontaneous dF/F')
pvals.SW.blank=permutationTest(SWtuned_NH.blank_rawDF,SWtuned_EN.blank_rawDF,10000);
pvals.CW.blank=permutationTest(CWtuned_NH.blank_rawDF,CWtuned_EN.blank_rawDF,10000);

% raw median dF/F responses
figure;
subplot(1,3,1);
hold on
plot_4cdfs(SWtuned_NH.BW_rawDF,SWtuned_EN.BW_rawDF,CWtuned_NH.BW_rawDF,CWtuned_EN.BW_rawDF)
ylabel('')
xlabel('median dF/F')
title('BW responses')
pvals.SW.BWraw=permutationTest(SWtuned_NH.BW_rawDF,SWtuned_EN.BW_rawDF,10000);
pvals.CW.BWraw=permutationTest(CWtuned_NH.BW_rawDF,CWtuned_EN.BW_rawDF,10000);


subplot(1,3,2);
hold on
plot_4cdfs(SWtuned_NH.nonBW_rawDF,SWtuned_EN.nonBW_rawDF,CWtuned_NH.nonBW_rawDF,CWtuned_EN.nonBW_rawDF)
ylabel('')
xlabel('median dF/F')
title('non-BW responses')
pvals.SW.nonBWraw=permutationTest(SWtuned_NH.nonBW_rawDF,SWtuned_EN.nonBW_rawDF,10000);
pvals.CW.nonBWraw=permutationTest(CWtuned_NH.nonBW_rawDF,CWtuned_EN.nonBW_rawDF,10000);

subplot(1,3,3);
hold on
plot_4cdfs(SWtuned_NH.CW_rawDF,SWtuned_EN.CW_rawDF,CWtuned_NH.CW_rawDF,CWtuned_EN.CW_rawDF)
legend('SW tuned NH','SW tuned EN','CW tuned NH','CW tuned EN')
ylabel('')
xlabel('median dF/F')
title('CW responses')
pvals.SW.CWraw=permutationTest(SWtuned_NH.CW_rawDF,SWtuned_EN.CW_rawDF,10000);
pvals.CW.CWraw=permutationTest(CWtuned_NH.CW_rawDF,CWtuned_EN.CW_rawDF,10000);

% Z-scored median dF/F responses
figure;
subplot(1,3,1);
hold on
plot_4cdfs(SWtuned_NH.BW_Z,SWtuned_EN.BW_Z,CWtuned_NH.BW_Z,CWtuned_EN.BW_Z)
ylabel('')
xlabel('Z-scored median')
title('BW responses')

subplot(1,3,2);
hold on
plot_4cdfs(SWtuned_NH.nonBW_Z,SWtuned_EN.nonBW_Z,CWtuned_NH.nonBW_Z,CWtuned_EN.nonBW_Z)
ylabel('')
xlabel('Z-scored median')
title('non-BW responses')

subplot(1,3,3);
hold on
plot_4cdfs(SWtuned_NH.CW_Z,SWtuned_EN.CW_Z,CWtuned_NH.CW_Z,CWtuned_EN.CW_Z)
legend('MW tuned NH','MW tuned EN','OW tuned NH','OW tuned EN')
ylabel('')
xlabel('Z-scored median')
title('CW responses')




    

end

