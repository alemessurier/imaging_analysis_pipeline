function [pvals_L4_tuning,pvals_L23_tuning,pvals_L4_point,pvals_L23_point]=somatotopy_plots_L4vL23()
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


% L4 - restrict ROIs to those in C1,C2,D1,D2
[ paths_EN,paths_NH ] = scnn1_paths;
whisk={'c1','c2','c3','d1','d2','d3','e1','e2','e3'};

[~,ROIs_NH,CWs_all_NH]=plot_allROI_positions( paths_NH,whisk,'um',1,2 );
[~,ROIs_EN,CWs_all_EN]=plot_allROI_positions( paths_EN,whisk,'um',1,2 );

for i=1:numel(ROIs_EN)
    inds_CRow=cellfun(@(x)ismember(x,{'c1','c2','c3'}),CWs_all_EN{i},'un',1);
    ROIs_CRow_EN{i}=ROIs_EN{i}(inds_CRow);
end

for i=1:numel(ROIs_NH)
    inds_CRow=cellfun(@(x)ismember(x,{'c1','c2','c3'}),CWs_all_NH{i},'un',1);
    ROIs_CRow_NH{i}=ROIs_NH{i}(inds_CRow);
end

    useEN=cellfun(@(x)length(x)>2,ROIs_CRow_EN);
    paths_EN_L4=paths_EN(useEN);
    ROIs_CRow_EN=ROIs_CRow_EN(useEN);


    useNH=cellfun(@(x)length(x)>2,ROIs_CRow_NH);
    paths_NH_L4=paths_NH(useNH);
    ROIs_CRow_NH=ROIs_CRow_NH(useNH);

    %%
 [tunedInCol_BW_NH,tunedInCol_BW_EN,percentMatch_BW_NH,percentMatch_BW_EN]=make_tunedInColPlots(paths_EN_L4,paths_NH_L4,1,2,ROIs_CRow_EN,ROIs_CRow_NH);
x=table([sum(tunedInCol_BW_NH); length(tunedInCol_BW_NH)-sum(tunedInCol_BW_NH)],...
    [sum(tunedInCol_BW_EN); length(tunedInCol_BW_EN)-sum(tunedInCol_BW_EN)],'VariableNames',...
    {'NH','EN'},'RowNames',{'CWtuned','notCWtuned'})
[h,p,stats] = fishertest(x)

meanPercMatch_BW_EN=mean(percentMatch_BW_EN)
semPercMatch_BW_EN=std(percentMatch_BW_EN)/sqrt(length(percentMatch_BW_EN))

meanPercMatch_BW_NH=mean(percentMatch_BW_NH)
semPercMatch_BW_NH=std(percentMatch_BW_NH)/sqrt(length(percentMatch_BW_NH))

%%
[ paths_EN,paths_NH ] = paths_to_include();
npSub=0;

[tunedInCol_EBW_NH,tunedInCol_EBW_EN,percentMatch_EBW_NH,percentMatch_EBW_EN]=make_tunedInColPlots(paths_EN,paths_NH,npSub,0,[],[]);
x=table([sum(tunedInCol_EBW_NH); length(tunedInCol_EBW_NH)-sum(tunedInCol_EBW_NH)],...
    [sum(tunedInCol_EBW_EN); length(tunedInCol_EBW_EN)-sum(tunedInCol_EBW_EN)],'VariableNames',...
    {'NH','EN'},'RowNames',{'CWtuned','notCWtuned'})
[h,p,stats] = fishertest(x);

meanPercMatch_EBW_EN=mean(percentMatch_EBW_EN)
semPercMatch_EBW_EN=std(percentMatch_EBW_EN)/sqrt(length(percentMatch_EBW_EN))

meanPercMatch_EBW_NH=mean(percentMatch_EBW_NH)
semPercMatch_EBW_NH=std(percentMatch_EBW_NH)/sqrt(length(percentMatch_EBW_NH))


%% plot percent of cells tuned to CW by radius (EBW)

% L23
[ paths_EN,paths_NH ] = paths_to_include();
npSub=0;
[mean_dist_NH_L23,mean_percMatch_NH_L23,sem_percMatch_NH_L23,percMatch_byBarr_NH_L23,percMatch_NH_L23,num_ROIs_NH_L23,CWtunedInds_binned_NH_L23,binEdges]=plot_CWtunedByRad(paths_NH,0,[],15,[]);
[mean_dist_EN_L23,mean_percMatch_EN_L23,sem_percMatch_EN_L23,percMatch_byBarr_EN_L23, percMatch_EN_L23,num_ROIs_EN_L23,CWtunedInds_binned_EN_L23]=plot_CWtunedByRad(paths_EN,0,[],binEdges,[]);


%L4

[mean_dist_NH_L4,mean_percMatch_NH_L4,sem_percMatch_NH_L4,percMatch_byBarr_NH_L4,percMatch_NH_L4,num_ROIs_NH_L4,CWtunedInds_binned_NH_L4,binEdges]=plot_CWtunedByRad(paths_NH_L4,1,2,binEdges,ROIs_CRow_NH);
[mean_dist_EN_L4,~,sem_percMatch_EN_L4,percMatch_byBarr_EN_L4, percMatch_EN_L4,num_ROIs_EN_L4,CWtunedInds_binned_EN_L4]=plot_CWtunedByRad(paths_EN_L4,1,2,binEdges,ROIs_CRow_EN);


P=zeros(1,15);%numel(CWtunedInds_binned_NH_L4));
for i=1:15%numel(CWtunedInds_binned_NH_L4)
    x=table([sum(CWtunedInds_binned_NH_L23{i}); length(CWtunedInds_binned_NH_L23{i})-sum(CWtunedInds_binned_NH_L23{i})],...
    [sum(CWtunedInds_binned_EN_L23{i}); length(CWtunedInds_binned_EN_L23{i})-sum(CWtunedInds_binned_EN_L23{i})],'VariableNames',...
    {'CT','EN'},'RowNames',{'CWtuned','notCWtuned'});
[h,P(i),stats] = fishertest(x);

end

pvals_L23_tuning.CWtunedByDist_EBW=P;

P=zeros(1,numel(CWtunedInds_binned_NH_L4));
for i=1:numel(CWtunedInds_binned_NH_L4)
    x=table([sum(CWtunedInds_binned_NH_L4{i}); length(CWtunedInds_binned_NH_L4{i})-sum(CWtunedInds_binned_NH_L4{i})],...
    [sum(CWtunedInds_binned_EN_L4{i}); length(CWtunedInds_binned_EN_L4{i})-sum(CWtunedInds_binned_EN_L4{i})],'VariableNames',...
    {'CT','EN'},'RowNames',{'CWtuned','notCWtuned'});
[h,P(i),stats] = fishertest(x);

end
pvals_L4_tuning.CWtunedByDist_EBW=P;


% compare first bin in L4 and L23 w/ fisher's test
for i=1:15
x=table([sum(CWtunedInds_binned_NH_L23{i}); length(CWtunedInds_binned_NH_L23{i})-sum(CWtunedInds_binned_NH_L23{i})],...
    [sum(CWtunedInds_binned_NH_L4{i}); length(CWtunedInds_binned_NH_L4{i})-sum(CWtunedInds_binned_NH_L4{i})],'VariableNames',...
    {'L23','L4'},'RowNames',{'CWtuned','notCWtuned'});
[h,p(i),stats] = fishertest(x);
end

P=zeros(1,numel(CWtunedInds_binned_NH_L23));
for i=1:numel(CWtunedInds_binned_NH_L23)
    P(i)=permutationTest(CWtunedInds_binned_NH_L23{i},CWtunedInds_binned_NH_L4{i},10000);
end
figure; hold on
plot(mean_dist_NH_L23,percMatch_NH_L23,'ko-','MarkerSize',4);
plot(mean_dist_EN_L23,percMatch_EN_L23,'ro-','MarkerSize',4);
plot(mean_dist_NH_L4,percMatch_NH_L4,'b*-');
plot(mean_dist_EN_L4,percMatch_EN_L4,'r*:');

xlabel('distance from CW center (um)')
ylabel('% tuned to CW (EBW)')

legend('NH L2/3','EN L2/3','NH L4','EN L4')

sem_NH=zeros(1,length(percMatch_NH_L23));
sem_EN=zeros(1,length(percMatch_EN_L23));
[~, ax] = totesComboPlot( mean_dist_NH_L23,percMatch_NH_L23,sem_NH,num_ROIs_NH_L23,[],mean_dist_EN_L23,percMatch_EN_L23,sem_EN,num_ROIs_EN_L23,[] );
ax(2).XLabel.String='distance from CW center (um)';
ax(2).YLabel.String='% tuned to CW (EBW)';
title('L23')
legend('NH','EN')

sem_NH=zeros(1,length(percMatch_NH_L4));
sem_EN=zeros(1,length(percMatch_EN_L4));
[~, ax] = totesComboPlot( mean_dist_NH_L4,percMatch_NH_L4,sem_NH,num_ROIs_NH_L4,[],mean_dist_EN_L4,percMatch_EN_L4,sem_EN,num_ROIs_EN_L4,[] );
ax(2).XLabel.String='distance from CW center (um)';
ax(2).YLabel.String='% tuned to CW (EBW)';
title('L4')
legend('NH','EN')


%% plot percent of cells tuned to CW by radius (absolute BW)

% L23
[ paths_EN,paths_NH ] = paths_to_include();
npSub=0;
[mean_dist_NH_L23,mean_percMatch_NH_L23,sem_percMatch_NH_L23,percMatch_byBarr_NH_L23,...
    percMatch_NH_L23,num_ROIs_NH_L23,CWtunedInds_binned_NH_L23,binEdges]=plot_CWtunedByRad_num(paths_NH,0,0,15,[]);
[mean_dist_EN_L23,mean_percMatch_EN_L23,sem_percMatch_EN_L23,percMatch_byBarr_EN_L23,...
    percMatch_EN_L23,num_ROIs_EN_L23,CWtunedInds_binned_EN_L23]=plot_CWtunedByRad_num(paths_EN,0,0,binEdges,[]);


% L4
% [ paths_EN,paths_NH ] = scnn1_paths;
npSub=1;
[mean_dist_NH_L4,mean_percMatch_NH_L4,sem_percMatch_NH_L4,percMatch_byBarr_NH_L4,...
    percMatch_NH_L4,num_ROIs_NH_L4,CWtunedInds_binned_NH_L4,binEdges]=plot_CWtunedByRad_num(paths_NH_L4,1,2,binEdges,ROIs_CRow_NH);
[mean_dist_EN_L4,mean_percMatch_EN_L4,sem_percMatch_EN_L4,percMatch_byBarr_EN_L4,...
    percMatch_EN_L4,num_ROIs_EN_L4,CWtunedInds_binned_EN_L4]=plot_CWtunedByRad_num(paths_EN_L4,1,2,binEdges,ROIs_CRow_EN);

P=zeros(1,10);
for i=1:10%numel(CWtunedInds_binned_NH_L4)
    P(i)=permutationTest(CWtunedInds_binned_NH_L4{i},CWtunedInds_binned_EN_L4{i},10000);
end

pvals_L4_tuning.CWtunedByDist_num=P;

P=zeros(1,numel(CWtunedInds_binned_NH_L23));
for i=1:numel(CWtunedInds_binned_NH_L23)
    P(i)=permutationTest(CWtunedInds_binned_NH_L23{i},CWtunedInds_binned_EN_L23{i},10000);
end

pvals_L23_tuning.CWtunedByDist_num=P;

P=zeros(1,10);
for i=1:10%numel(CWtunedInds_binned_NH_L4)
    P(i)=permutationTest(CWtunedInds_binned_NH_L4{i},CWtunedInds_binned_NH_L23{i},10000);
end

pvals_L23vL4_tuning.CWtunedByDist_num=P;

figure; hold on
plot(mean_dist_NH_L23,percMatch_NH_L23,'ko-','MarkerSize',4);
plot(mean_dist_EN_L23,percMatch_EN_L23,'ro-','MarkerSize',4);
plot(mean_dist_NH_L4,percMatch_NH_L4,'b*-');
plot(mean_dist_EN_L4,percMatch_EN_L4,'r*:');

xlabel('distance from CW center (um)')
ylabel('% tuned to CW (numerical)')
legend('NH L2/3','EN L2/3','NH L4','EN L4')

sem_NH=zeros(1,length(percMatch_NH_L23));
sem_EN=zeros(1,length(percMatch_EN_L23));
[~, ax] = totesComboPlot( mean_dist_NH_L23,percMatch_NH_L23,sem_NH,num_ROIs_NH_L23,[],mean_dist_EN_L23,percMatch_EN_L23,sem_EN,num_ROIs_EN_L23,[] );
ax(2).XLabel.String='distance from CW center (um)';
ax(2).YLabel.String='% tuned to CW (numerical)';
title('L23')
legend('NH','EN')

sem_NH=zeros(1,length(percMatch_NH_L4));
sem_EN=zeros(1,length(percMatch_EN_L4));
[~, ax] = totesComboPlot( mean_dist_NH_L4,percMatch_NH_L4,sem_NH,num_ROIs_NH_L4,[],mean_dist_EN_L4,percMatch_EN_L4,sem_EN,num_ROIs_EN_L4,[] );
ax(2).XLabel.String='distance from CW center (um)';
ax(2).YLabel.String='% tuned to CW (numerical)';
title('L4')
legend('NH','EN')

P=zeros(1,15);%numel(CWtunedInds_binned_NH_L4));
for i=1:15%numel(CWtunedInds_binned_NH_L4)
    x=table([sum(CWtunedInds_binned_NH_L23{i}); length(CWtunedInds_binned_NH_L23{i})-sum(CWtunedInds_binned_NH_L23{i})],...
    [sum(CWtunedInds_binned_EN_L23{i}); length(CWtunedInds_binned_EN_L23{i})-sum(CWtunedInds_binned_EN_L23{i})],'VariableNames',...
    {'CT','EN'},'RowNames',{'CWtuned','notCWtuned'});
[h,P(i),stats] = fishertest(x);

end

pvals_L23_tuning.CWtunedByDist_num=P;

P=zeros(1,numel(CWtunedInds_binned_NH_L4));
for i=1:numel(CWtunedInds_binned_NH_L4)
    x=table([sum(CWtunedInds_binned_NH_L4{i}); length(CWtunedInds_binned_NH_L4{i})-sum(CWtunedInds_binned_NH_L4{i})],...
    [sum(CWtunedInds_binned_EN_L4{i}); length(CWtunedInds_binned_EN_L4{i})-sum(CWtunedInds_binned_EN_L4{i})],'VariableNames',...
    {'CT','EN'},'RowNames',{'CWtuned','notCWtuned'});
[h,P(i),stats] = fishertest(x);

end
pvals_L4_tuning.CWtunedByDist_num=P;

% compare bins in L4 and L23 w/ fisher's test
p=zeros(15,1);
for i=1:15
x=table([sum(CWtunedInds_binned_NH_L23{i}); length(CWtunedInds_binned_NH_L23{i})-sum(CWtunedInds_binned_NH_L23{i})],...
    [sum(CWtunedInds_binned_NH_L4{i}); length(CWtunedInds_binned_NH_L4{i})-sum(CWtunedInds_binned_NH_L4{i})],'VariableNames',...
    {'L23','L4'},'RowNames',{'CWtuned','notCWtuned'});
[h,p(i),stats] = fishertest(x);
end

%%
% L2/3
[ paths_EN,paths_NH ] = paths_to_include();
[ax_L23,pvals_L23_point,CWresp_byDist_L23]=plotPointRep(paths_EN,paths_NH,0,15,[],[]);


% L4

% [ paths_EN,paths_NH ] = scnn1_paths;
[ax_L4,pvals_L4_point,CWresp_byDist_L4]=plotPointRep(paths_EN_L4,paths_NH_L4,1,15,ROIs_CRow_EN,ROIs_CRow_NH);

boundsY=ax_L4(2).YLim;
boundsX=ax_L23(2).XLim;

ax_L4(2).Title.String='Layer 4';
ax_L4(2).XLim=boundsX;
ax_L23(2).YLim=boundsY;
ax_L23(2).Title.String='L23';


%% RF size
[ paths_EN,paths_NH ] = paths_to_include();
compareRFsize(paths_EN,paths_NH,0,'L23',[],[])

% [ paths_EN,paths_NH ] = scnn1_paths;
compareRFsize(paths_EN_L4,paths_NH_L4,1,'L4',ROIs_CRow_EN,ROIs_CRow_NH)

%% find percent responsive

[ percent_responsive_L23_NH,resp_byField_L23_NH ] = find_percent_responsive( paths_NH,0 )
[ percent_responsive_L23_EN,resp_byField_L23_EN ] = find_percent_responsive( paths_EN,0 )

[ percent_responsive_L4_NH,resp_byField_L4_NH ] = find_percent_responsive( paths_NH_L4,1 )
[ percent_responsive_L4_EN,resp_byField_L4_EN ] = find_percent_responsive( paths_EN_L4,1 )

%% compare RF sharpness
[ paths_EN,paths_NH ] = paths_to_include();
p_L23_BW=compareRFsharpness(paths_EN,paths_NH,0,'L23',[],[],'BW')

% [ paths_EN,paths_NH ] = scnn1_paths;
p_L4_BW=compareRFsharpness(paths_EN_L4,paths_NH_L4,1,'L4',ROIs_CRow_EN,ROIs_CRow_NH,'BW')

% [ paths_EN,paths_NH ] = paths_to_include();
p_L23_CW=compareRFsharpness(paths_EN,paths_NH,0,'L23',[],[],'CW')

% [ paths_EN,paths_NH ] = scnn1_paths;
p_L4_CW=compareRFsharpness(paths_EN_L4,paths_NH_L4,1,'L4',ROIs_CRow_EN,ROIs_CRow_NH,'CW')

%%
    
save(['E:\figures20180528\figure1_somatotopy\L23vL4\pvals_',date,'.mat'],'pvals_L23_point','pvals_L23_tuning','pvals_L4_point','pvals_L4_tuning')

   

    


   





    
end