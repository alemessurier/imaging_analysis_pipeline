function voronoi_decode_detect(paths_EN,paths_NH)

ROI_pos_plot_EN=cell(length(paths_EN),1);
perfByWhisk_EN=ROI_pos_plot_EN;
PbyWshuff_EN=ROI_pos_plot_EN;

ROI_pos_plot_NH=cell(length(paths_NH),1);
perfByWhisk_NH=ROI_pos_plot_NH;
PbyWshuff_NH=ROI_pos_plot_NH;

      for K=1:length(paths_EN)
           [ ROI_pos_plot_EN{K},perfByWhisk_EN{K},PbyWshuff_EN{K} ] = decoder_PandPos( paths_EN{K},'all' );
        end
        
      for K=1:length(paths_NH)
           [ ROI_pos_plot_NH{K},perfByWhisk_NH{K},PbyWshuff_NH{K} ] = decoder_PandPos( paths_NH{K},'all' );
        end

ROI_pos_plot_EN=cat(1,ROI_pos_plot_EN{:});
perfByWhisk_EN=cat(1,perfByWhisk_EN{:});
PbyWshuff_EN=cat(1,PbyWshuff_EN{:});

ROI_pos_plot_NH=cat(1,ROI_pos_plot_NH{:});
perfByWhisk_NH=cat(1,perfByWhisk_NH{:});
PbyWshuff_NH=cat(1,PbyWshuff_NH{:});

%% kmeans cluster points, voronoi diagram - response size vs. um distance

pool = gcp;                      % Invokes workers  
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);


nClust=500;
numBins=15;

% ROI_pos_plot_NH=fliplr(ROI_pos_plot_NH);
% ROI_pos_plot_EN=fliplr(ROI_pos_plot_EN);

%find cmap bounds
[pos_idx,ROIpos_clust_NH]=kmeans(ROI_pos_plot_NH,nClust,'Options',options,'MaxIter',10000,...
    'Display','final','Replicates',10);
for i=1:nClust
    Rs=perfByWhisk_NH(pos_idx==i);
    cVar_clust_NH(i)=mean(Rs);
    Rshuff=PbyWshuff_NH(pos_idx==i);
    cShuff_clust_NH(i)=mean(Rshuff);
end

[pos_idx,ROIpos_clust_EN]=kmeans(ROI_pos_plot_EN,nClust,'Options',options,'MaxIter',10000,...
    'Display','final','Replicates',10);
for i=1:nClust
     Rs=perfByWhisk_EN(pos_idx==i);
    cVar_clust_EN(i)=mean(Rs);
    Rshuff=PbyWshuff_EN(pos_idx==i);
    cShuff_clust_EN(i)=mean(Rshuff);
end

% cluster, make voronoi plots for non-shuffled performance by ROI
cmapBounds(1)=min([cVar_clust_EN,cVar_clust_NH]);
cmapBounds(2)=max([cVar_clust_EN,cVar_clust_NH]);

figure; hold on
subplot(1,2,1)
voronoi_ROI_pos(ROIpos_clust_NH,cVar_clust_NH,'numBins',numBins,'cmapBounds',cmapBounds);
tmp=viscircles([0 0],150,'EdgeColor','k','LineStyle','--','LineWidth',2);
tmp.Children(1).Marker='none';
title('NH')

subplot(1,2,2)
voronoi_ROI_pos(ROIpos_clust_EN,cVar_clust_EN,'numBins',numBins,'cmapBounds',cmapBounds);
tmp2=viscircles([0 0],150,'EdgeColor','k','LineStyle','--','LineWidth',2);
tmp2.Children(1).Marker='none';
title('EN')

% shuffled performance
cmapBounds(1)=min([cShuff_clust_EN,cShuff_clust_NH]);
cmapBounds(2)=max([cShuff_clust_EN,cShuff_clust_NH]);

figure; hold on
subplot(1,2,1)
voronoi_ROI_pos(ROIpos_clust_NH,cShuff_clust_NH,'numBins',numBins,'cmapBounds',cmapBounds);
tmp=viscircles([0 0],150,'EdgeColor','k','LineStyle','--','LineWidth',2);
tmp.Children(1).Marker='none';
title('NH')

subplot(1,2,2)
voronoi_ROI_pos(ROIpos_clust_EN,cShuff_clust_EN,'numBins',numBins,'cmapBounds',cmapBounds);
tmp2=viscircles([0 0],150,'EdgeColor','k','LineStyle','--','LineWidth',2);
tmp2.Children(1).Marker='none';
title('EN')

end