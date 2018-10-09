function compare_npCorrection_20171327()

% %% plot fraction of ROIs reponsive to whisker deflection for different conditions
% 
 [ paths_EN,paths_NH ] = paths_to_include();
%%
for K=1:length(paths_NH)
    [~,traceByStim,~,~,permTestResults] = load_nonNPsub_data( paths_NH{K} );
    [~,~,sigInds{K}]=find_sigROIs(permTestResults,traceByStim);
end

sigNH_noNP_byField=cellfun(@(x)sum(x)/length(x),sigInds);
sigInds=cat(2,sigInds{:});
sig_NH_noNPsub=sum(sigInds)/length(sigInds);

clear sigInds

for K=1:length(paths_EN)
    [~,traceByStim,~,~,permTestResults] = load_nonNPsub_data( paths_EN{K} );
    [~,~,sigInds{K}]=find_sigROIs(permTestResults,traceByStim);
end
sigEN_noNP_byField=cellfun(@(x)sum(x)/length(x),sigInds)
sigInds=cat(2,sigInds{:});
sig_EN_noNPsub=sum(sigInds)/length(sigInds);

clear sigInds

for K=1:length(paths_NH)
    [~,traceByStim,~,~,permTestResults] = load_NPsub_data( paths_NH{K} );
    [~,~,sigInds{K}]=find_sigROIs(permTestResults,traceByStim);
end
sigNH_NPsub_byField=cellfun(@(x)sum(x)/length(x),sigInds)
sigInds=cat(2,sigInds{:});
sig_NH_NPsub=sum(sigInds)/length(sigInds);

clear sigInds

for K=1:length(paths_EN)
    [~,traceByStim,~,~,permTestResults] = load_NPsub_data( paths_EN{K} );
    [~,~,sigInds{K}]=find_sigROIs(permTestResults,traceByStim);
end

sigEN_NPsub_byField=cellfun(@(x)sum(x)/length(x),sigInds)
sigInds=cat(2,sigInds{:});
sig_EN_NPsub=sum(sigInds)/length(sigInds);

figure; hold on
H=bar([sig_NH_noNPsub,sig_EN_noNPsub,sig_NH_NPsub,sig_EN_NPsub]);
H.FaceColor='none';
H.LineWidth=1.5;
ax=gca;
ax.XTickLabel={'','NH','EN','NH,NP corrected','EN,NP corrected',''}
ax.XTickLabelRotation=45;
ylabel('fraction whisker responsive ROIs')

figure; hold on

pl(1)=notBoxPlot(sigNH_noNP_byField,1,[],'line')
pl(1).sd.Color='none';
pl(1).data.MarkerFaceColor='none';


pl(2)=notBoxPlot(sigEN_noNP_byField,2,[],'line')
pl(2).sd.Color='none';
pl(2).data.MarkerFaceColor='none';

pl(3)=notBoxPlot(sigNH_NPsub_byField,3,[],'line')
pl(3).sd.Color='none';
pl(3).data.MarkerFaceColor='none';

pl(4)=notBoxPlot(sigEN_NPsub_byField,4,[],'line')
pl(4).sd.Color='none';
pl(4).data.MarkerFaceColor='none';

ax2=gca;
ax2.XTick=1:4;
ax2.XLim=[0 5]
ax2.XTickLabels={'NH','EN','NH,NP corrected','EN,NP corrected'}
ax2.XTickLabelRotation=45;
ylabel('fraction of responsive ROIs/field')


%%
paths_all=[paths_NH,paths_EN];

for K=1:length(paths_all)
    load(strcat(paths_all{K},'whiskPref.mat'),'whiskPref');
    cellNames=fieldnames(whiskPref);
    numBW{K}=cellfun(@(x)length(whiskPref.(x){1}),cellNames,'Uni',1);
    numSigW{K}=cellfun(@(x)length(horzcat(whiskPref.(x){:})),cellNames,'Uni',1);
end

numBW_all=vertcat(numBW{:});
numSigW_all=vertcat(numSigW{:});

for K=1:length(paths_all)
    load(strcat(paths_all{K},'whiskPref_NPsub.mat'),'whiskPref');
    cellNames=fieldnames(whiskPref);
    numBW_NP{K}=cellfun(@(x)length(whiskPref.(x){1}),cellNames,'Uni',1);
    numSigW_NP{K}=cellfun(@(x)length(horzcat(whiskPref.(x){:})),cellNames,'Uni',1);
end
numBW_NP_all=vertcat(numBW_NP{:});
numSigW_NP_all=vertcat(numSigW_NP{:});

figure; hold on
cdfplot(numBW_all)
cdfplot(numBW_NP_all)
legend('no corr','corr')
ylabel('')
xlabel('# equivalent BWs')
[~,pval]=ttest2(numBW_all,numBW_NP_all);
title(strcat('p=',num2str(pval)))

figure; hold on
histogram(numBW_NP_all,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(numBW_all,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
ax=gca;
ax.XTick=1:9;
ax.XTickLabel=1:9;
% ax.Children.FaceColor='none';
% ax.Children(1).LineWidth=1.5;
% ax.Children(2).LineWidth=1.5;
legend('np corr','no corr')
xlabel('# equivalent BWs')
ylabel('fraction of ROIs')


figure; hold on
histogram(numSigW_NP_all,1:10,'Normalization','probability','FaceColor','k','FaceAlpha',0.6,'EdgeColor','k');%,'FaceAlpha',0.5);
histogram(numSigW_all,1:10,'Normalization','probability','FaceColor','r','FaceAlpha',0.4,'EdgeColor','r');
ax=gca;
ax.XTick=1:9;
ax.XTickLabel=1:9;
% ax.Children.FaceColor='none';
% ax.Children(1).LineWidth=1.5;
% ax.Children(2).LineWidth=1.5;
legend('np corr','no corr')
xlabel('# whiskers in RF')
ylabel('fraction of ROIs')
[~,p]=ttest2(numSigW_NP_all,numSigW_all)
title(['p=',num2str(p)]);
%% Find distributions of deltaF/F for different conditions

[meanDF_NH,stdDF_NH,skewDF_NH]=find_deltaF_stats(paths_NH,0)
[meanDF_EN,stdDF_EN,skewDF_EN]=find_deltaF_stats(paths_EN,0)
[meanDF_NH_NP,stdDF_NH_NP,skewDF_NH_NP]=find_deltaF_stats(paths_NH,1)
[meanDF_EN_NP,stdDF_EN_NP,skewDF_EN_NP]=find_deltaF_stats(paths_EN,1)


% find skewness, mean, stdev of dF/F for each ROI
means=figure; hold on
m(1)=cdfplot(meanDF_NH)
m(2)=cdfplot(meanDF_EN)
m(3)=cdfplot(meanDF_NH_NP)
m(4)=cdfplot(meanDF_EN_NP)
legend('NH','EN','NH,corr','EN,corr')
xlabel('mean dF/F (by ROI)')
ylabel('')
title('')

stds=figure; hold on
s(1)=cdfplot(stdDF_NH)
s(2)=cdfplot(stdDF_EN)
s(3)=cdfplot(stdDF_NH_NP)
s(4)=cdfplot(stdDF_EN_NP)
legend('NH','EN','NH,corr','EN,corr')
xlabel('stdev of dF/F (by ROI)')
ylabel('')
title('')


skew=figure; hold on
sk(1)=cdfplot(skewDF_NH)
sk(2)=cdfplot(skewDF_EN)
sk(3)=cdfplot(skewDF_NH_NP)
sk(4)=cdfplot(skewDF_EN_NP)
legend('NH','EN','NH,corr','EN,corr')
xlabel('skewness of dF/F (by ROI)')
ylabel('')
title('')
% stds=figure; hold on
% skew=figure; hold on



    function [meanDF_all,stdDF_all,skewDF_all]=find_deltaF_stats(paths,npSub)
        for K=1:length(paths)
            cd(paths{K});
            switch npSub
                case 0
                    fname_s1=dir('step1_*');
                case 1
                    fname_s1=dir('step1NP_*');
            end
            load(strcat(paths{K},fname_s1(end).name),'deltaF');
            fns=fieldnames(deltaF);
            cellNames=fieldnames(deltaF.(fns{1}));
            meanDF=zeros(1,length(cellNames));
            stdDF=zeros(1,length(cellNames));
            skewDF=zeros(1,length(cellNames));
            for i=1:length(cellNames)
                dFall=cellfun(@(x)deltaF.(x).(cellNames{i})(20:end-20),fns,'Uni',0);
                dFall=cat(1,dFall{:});
                meanDF(i)=mean(dFall);
                stdDF(i)=std(dFall);
                skewDF(i)=max(dFall)/norm(dFall);
            end
            %     figure(means)
            %     cdfplot(meanDF);
            %     figure(stds)
            %     cdfplot(stdDF)
            %     figure(skew)
            %     cdfplot(skewDF)
            meanDF_all{K}=meanDF;
            stdDF_all{K}=stdDF;
            skewDF_all{K}=skewDF;
            
            
        end
        meanDF_all=cat(2,meanDF_all{:});
        stdDF_all=cat(2,stdDF_all{:});
        skewDF_all=cat(2,skewDF_all{:});
    end
end
