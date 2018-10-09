function make_responseDistPlots_npSubCompare( paths_EN,paths_NH )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
for r=1:4
    [ blank_rawDF_EN,BW_rawDF_EN,nonBW_rawDF_EN,CW_rawDF_EN,BW_Z_EN,nonBW_Z_EN,CW_Z_EN,numSigW_EN]=find_response_dists(paths_EN,1,r,[]);
    [ blank_rawDF_NH,BW_rawDF_NH,nonBW_rawDF_NH,CW_rawDF_NH,BW_Z_NH,nonBW_Z_NH,CW_Z_NH,numSigW_NH]=find_response_dists(paths_NH,1,r,[]);
    
    
    
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
    title(strcat('spontaneous dF/F, p=',num2str(pval),', r=',num2str(r)))
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
    pval=permutationTest(BW_rawDF_EN,BW_rawDF_NH,100000);
    title(strcat('BW response, p=',num2str(pval),', r=',num2str(r)))
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
    title(strcat('non-BW response, p=',num2str(pval),', r=',num2str(r)))
    
    
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
    title(strcat('BW response, p=',num2str(pval),', r=',num2str(r)))
    % tmp=gca;
    % tmp.XLim(2)=0.15
    
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
    title(strcat('non-BW response, p=',num2str(pval),', r=',num2str(r)))
    
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
    title(strcat('CW response,p=',num2str(pval),', r=',num2str(r)))
    
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
    title(strcat('CW response,p=',num2str(pval),', r=',num2str(r)))
end
%% median response analysis for Enriched animals

% 
%     function [ blank_rawDF,BW_rawDF,nonBW_rawDF,CW_rawDF,BW_Z,nonBW_Z,CW_Z,numSigW]=find_response_dists(pathNames,r)
%         
%         
%         
%         
%         for K=1:length(pathNames)
%             
%                     [~,traceByStim,sponTrace,framesEvoked,permTestResults,~,ROIsInBarrel,~,~,Stimuli,deltaF,sampRate ] = load_NPsub_data( pathNames{K} );
%             
%                     deltaF=deltaF(r);
%                     traceByStim=traceByStim(r);
%                     sponTrace=sponTrace(r);
%                     permTestResults=permTestResults(r);
%             [ sponTraceRaw ] = make_sponTrace_noBLS( Stimuli,sampRate(1),deltaF,0.5,3 );
%             [sigCells,numSig]=find_sigROIs( permTestResults,traceByStim ); %only include significantly responsive cells
%             cellNames=fieldnames(traceByStim);
%             [ rawDF(K),Zscored(K) ] = activity_medians(sigCells,traceByStim,sponTraceRaw,sponTrace,framesEvoked );
%             
%             [ ~,Zscored_CW{K},medians_raw_CW{K} ] = find_CW_responses_Z(sigCells,ROIsInBarrel,traceByStim,sponTrace,framesEvoked);
%       
%             numSigW{K}=cellfun(@(x)sum(x),numSig,'Uni',1);
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
%       
%         numSigW=cat(2,numSigW{:});
%     end






end

