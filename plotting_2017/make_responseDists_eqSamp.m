%%


[ROIsByField_EN, ROIsByField_NH]=resample_equalPos( paths_NH,paths_EN,type,ref_barrel,reps,npSub );
%%
parfor i=1:reps
    [ blank_rawDF_EN{i},BW_rawDF_EN{i},nonBW_rawDF_EN{i},CW_rawDF_EN{i},...
        BW_Z_EN{i},nonBW_Z_EN{i},CW_Z_EN{i},...
        fano_by_whisk_EN{i},medianRs_EN{i},fano_CW_EN{i}]=find_response_dists(paths_EN,npSub,ROIsByField_EN{i});
   
    [ blank_rawDF_NH{i},BW_rawDF_NH{i},nonBW_rawDF_NH{i},CW_rawDF_NH{i},...
        BW_Z_NH{i},nonBW_Z_NH{i},CW_Z_NH{i},...
        fano_by_whisk_NH{i},medianRs_NH{i},fano_CW_NH{i}]=find_response_dists(paths_NH,npSub,ROIsByField_NH{i});
end

%% histograms of differences in means: blank

mean_spon_NH=cellfun(@mean,blank_rawDF_NH);
mean_spon_EN=cellfun(@mean,blank_rawDF_EN);
mean_spon_diff=mean_spon_NH-mean_spon_EN;

figure; hold on
hist(mean_spon_diff)
title('difference in mean blank response, NH-EN')

%% histograms of differences in means: CW

mean_CW_NH=cellfun(@mean,CW_rawDF_NH);
mean_CW_EN=cellfun(@mean,CW_rawDF_EN);
mean_CW_diff=mean_CW_NH-mean_CW_EN;

figure; hold on
hist(mean_CW_diff)
title('difference in mean CW response, NH-EN')

mean_CW_NH=cellfun(@mean,CW_Z_NH);
mean_CW_EN=cellfun(@mean,CW_Z_EN);
mean_CW_diff=mean_CW_NH-mean_CW_EN;

figure; hold on
hist(mean_CW_diff)
title('difference in mean Z-scored CW response, NH-EN')

%% histograms of differences in means: BW

mean_BW_NH=cellfun(@mean,BW_rawDF_NH);
mean_BW_EN=cellfun(@mean,BW_rawDF_EN);
mean_BW_diff=mean_BW_NH-mean_BW_EN;

figure; hold on
hist(mean_BW_diff)
title('difference in mean BW response, NH-EN')

mean_BW_NH=cellfun(@mean,BW_Z_NH);
mean_BW_EN=cellfun(@mean,BW_Z_EN);
mean_BW_diff=mean_BW_NH-mean_BW_EN;

figure; hold on
hist(mean_BW_diff)
title('difference in mean Z-scored BW response, NH-EN')

%% histograms of differences in means: BW

mean_nonBW_NH=cellfun(@mean,nonBW_rawDF_NH);
mean_nonBW_EN=cellfun(@mean,nonBW_rawDF_EN);
mean_nonBW_diff=mean_nonBW_NH-mean_nonBW_EN;

figure; hold on
hist(mean_nonBW_diff)
title('difference in mean nonBW response, NH-EN')

mean_nonBW_NH=cellfun(@mean,nonBW_Z_NH);
mean_nonBW_EN=cellfun(@mean,nonBW_Z_EN);
mean_nonBW_diff=mean_nonBW_NH-mean_nonBW_EN;

figure; hold on
hist(mean_nonBW_diff)
title('difference in mean Z-scored nonBW response, NH-EN')
%% histograms of differences in means: fano to CW

fano_CW_NH=cellfun(@(x)cat(2,x{:}),fano_CW_NH,'un',0);
fano_CW_EN=cellfun(@(x)cat(2,x{:}),fano_CW_EN,'un',0);

mean_fano_CW_NH=cellfun(@mean,fano_CW_NH);
mean_fano_CW_EN=cellfun(@mean,fano_CW_EN);
mean_fano_CW_diff=mean_fano_CW_NH-mean_fano_CW_EN;


figure; hold on
hist(mean_fano_CW_diff)
title('difference in mean ff to CW, NH-EN')

figure; hold on
plot_2cdfs(mean_fano_CW_NH,mean_fano_CW_EN)

%% manually sub-sample: only cells in C1,C2

[~,ROIs_NH,CWs_all_NH]=plot_allROI_positions( paths_NH,whisk,type,npSub );
[~,ROIs_EN,CWs_all_EN]=plot_allROI_positions( paths_EN,whisk,type,npSub );

for i=1:numel(ROIs_EN)
    inds_CRow=cellfun(@(x)ismember(x,{'c1','c2'}),CWs_all_EN{i},'un',1);
    ROIs_CRow_EN{i}=ROIs_EN{i}(inds_CRow);
end

for i=1:numel(ROIs_NH)
    inds_CRow=cellfun(@(x)ismember(x,{'c1','c2'}),CWs_all_NH{i},'un',1);
    ROIs_CRow_NH{i}=ROIs_NH{i}(inds_CRow);
end

make_responseDistPlots_eqSamp( paths_EN,paths_NH,npSub,ROIs_CRow_EN,ROIs_CRow_NH )