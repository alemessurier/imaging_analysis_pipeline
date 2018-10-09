function [pvals_L23_tuning,pvals_L23_point]=somatotopy_plots_L23_npSubCompare_2018(rval,ROIs_NH,ROIs_EN)
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


% %%
% [ paths_EN,paths_NH ] = paths_to_include();
% 
% 
% [tunedInCol_EBW_NH,tunedInCol_EBW_EN,percentMatch_EBW_NH,percentMatch_EBW_EN]=make_tunedInColPlots(paths_EN,paths_NH,npSub,[],[]);
% x=table([sum(tunedInCol_EBW_NH); length(tunedInCol_EBW_NH)-sum(tunedInCol_EBW_NH)],...
%     [sum(tunedInCol_EBW_EN); length(tunedInCol_EBW_EN)-sum(tunedInCol_EBW_EN)],'VariableNames',...
%     {'NH','EN'},'RowNames',{'CWtuned','notCWtuned'})
% [h,p,stats] = fishertest(x);
% 
% meanPercMatch_EBW_EN=mean(percentMatch_EBW_EN)
% semPercMatch_EBW_EN=std(percentMatch_EBW_EN)/sqrt(length(percentMatch_EBW_EN))
% 
% meanPercMatch_EBW_NH=mean(percentMatch_EBW_NH)
% semPercMatch_EBW_NH=std(percentMatch_EBW_NH)/sqrt(length(percentMatch_EBW_NH))
% 
% 

% L23
% [ paths_EN,paths_NH ] = paths_to_include();
% for K=1:length(paths_NH)
%     [traceByStim,~,framesEvoked,permTestResults] = load_NPsub_data_L23( paths_NH{K},1 )
%     if isempty(ROIs_NH)
%         ROIs_NH{K}=find_sigROIs(permTestResults,traceByStim);
%     end
% end
% 
% for K=1:length(paths_EN)
%     [traceByStim,~,framesEvoked,permTestResults] = load_NPsub_data_L23( paths_EN{K},1 )
%     if isempty(ROIs_EN)
%         ROIs_EN{K}=find_sigROIs(permTestResults,traceByStim);
%     end
% end
%% plot percent of cells tuned to CW by radius (EBW)

[ paths_EN,paths_NH ] = paths_to_include();

npSub=0;
[mean_dist_NH_L23,mean_percMatch_NH_L23,sem_percMatch_NH_L23,percMatch_byBarr_NH_L23,percMatch_NH_L23,num_ROIs_NH_L23,CWtunedInds_binned_NH_L23,binEdges]=plot_CWtunedByRad(paths_NH,1,rval,15,ROIs_NH);
[mean_dist_EN_L23,mean_percMatch_EN_L23,sem_percMatch_EN_L23,percMatch_byBarr_EN_L23, percMatch_EN_L23,num_ROIs_EN_L23,CWtunedInds_binned_EN_L23]=plot_CWtunedByRad(paths_EN,1,rval,binEdges,ROIs_EN);




P=zeros(1,15);%numel(CWtunedInds_binned_NH_L4));
for i=1:15%numel(CWtunedInds_binned_NH_L4)
    x=table([sum(CWtunedInds_binned_NH_L23{i}); length(CWtunedInds_binned_NH_L23{i})-sum(CWtunedInds_binned_NH_L23{i})],...
    [sum(CWtunedInds_binned_EN_L23{i}); length(CWtunedInds_binned_EN_L23{i})-sum(CWtunedInds_binned_EN_L23{i})],'VariableNames',...
    {'CT','EN'},'RowNames',{'CWtuned','notCWtuned'});
[h,P(i),stats] = fishertest(x);

end

pvals_L23_tuning.CWtunedByDist_EBW=P;



sem_NH=zeros(1,length(percMatch_NH_L23));
sem_EN=zeros(1,length(percMatch_EN_L23));
[~, ax] = totesComboPlot( mean_dist_NH_L23,percMatch_NH_L23,sem_NH,num_ROIs_NH_L23,[],mean_dist_EN_L23,percMatch_EN_L23,sem_EN,num_ROIs_EN_L23,[] );
ax(2).XLabel.String='distance from CW center (um)';
ax(2).YLabel.String='% tuned to CW (EBW)';
title('L23')
legend('NH','EN')



%% plot percent of cells tuned to CW by radius (absolute BW)

% L23
[ paths_EN,paths_NH ] = paths_to_include();
% paths_NH=paths_NH(1:end-1);
npSub=0;
[mean_dist_NH_L23,mean_percMatch_NH_L23,sem_percMatch_NH_L23,percMatch_byBarr_NH_L23,percMatch_NH_L23,num_ROIs_NH_L23,CWtunedInds_binned_NH_L23,binEdges]=plot_CWtunedByRad_num(paths_NH,1,rval,15,ROIs_NH);
[mean_dist_EN_L23,mean_percMatch_EN_L23,sem_percMatch_EN_L23,percMatch_byBarr_EN_L23, percMatch_EN_L23,num_ROIs_EN_L23,CWtunedInds_binned_EN_L23]=plot_CWtunedByRad_num(paths_EN,1,rval,binEdges,ROIs_EN);


P=zeros(1,numel(CWtunedInds_binned_NH_L23));
for i=1:numel(CWtunedInds_binned_NH_L23)
    P(i)=permutationTest(CWtunedInds_binned_NH_L23{i},CWtunedInds_binned_EN_L23{i},10000);
end

pvals_L23_tuning.CWtunedByDist_num=P;


sem_NH=zeros(1,length(percMatch_NH_L23));
sem_EN=zeros(1,length(percMatch_EN_L23));
[~, ax] = totesComboPlot( mean_dist_NH_L23,percMatch_NH_L23,sem_NH,num_ROIs_NH_L23,[],mean_dist_EN_L23,percMatch_EN_L23,sem_EN,num_ROIs_EN_L23,[] );
ax(2).XLabel.String='distance from CW center (um)';
ax(2).YLabel.String='% tuned to CW (numerical)';
title('L23')
legend('NH','EN')


%%
% L2/3
[ paths_EN,paths_NH ] = paths_to_include();
% paths_NH=paths_NH(1:end-1);
[ax_L23,pvals_L23_point,CWresp_byDist_L23]=plotPointRep(paths_EN,paths_NH,1,rval,15,ROIs_EN,ROIs_NH);




% %% RF size
% [ paths_EN,paths_NH ] = paths_to_include();
% compareRFsize(paths_EN,paths_NH,0,'L23',[],[])
% 
% % [ paths_EN,paths_NH ] = scnn1_paths;
% compareRFsize(paths_EN_L4,paths_NH_L4,1,'L4',ROIs_CRow_EN,ROIs_CRow_NH)
% 
% %% find percent responsive
% 
% [ percent_responsive_L23_NH,resp_byField_L23_NH ] = find_percent_responsive( paths_NH,0 )
% [ percent_responsive_L23_EN,resp_byField_L23_EN ] = find_percent_responsive( paths_EN,0 )
% 
% [ percent_responsive_L4_NH,resp_byField_L4_NH ] = find_percent_responsive( paths_NH_L4,1 )
% [ percent_responsive_L4_EN,resp_byField_L4_EN ] = find_percent_responsive( paths_EN_L4,1 )
% 
% %% compare RF sharpness
% [ paths_EN,paths_NH ] = paths_to_include();
% p_L23_BW=compareRFsharpness(paths_EN,paths_NH,0,'L23',[],[],'BW')
% 
% [ paths_EN,paths_NH ] = scnn1_paths;
% p_L4_BW=compareRFsharpness(paths_EN_L4,paths_NH_L4,1,'L4',ROIs_CRow_EN,ROIs_CRow_NH,'BW')
% 
% [ paths_EN,paths_NH ] = paths_to_include();
% p_L23_CW=compareRFsharpness(paths_EN,paths_NH,0,'L23',[],[],'CW')
% 
% [ paths_EN,paths_NH ] = scnn1_paths;
% p_L4_CW=compareRFsharpness(paths_EN_L4,paths_NH_L4,1,'L4',ROIs_CRow_EN,ROIs_CRow_NH,'CW')
% 
% %%
%     
% save(['E:\figures20180528\figure1_somatotopy\L23vL4\pvals_',date,'.mat'],'pvals_L23_point','pvals_L23_tuning','pvals_L4_point','pvals_L4_tuning')

   

    


   





    
end