function voronoiPlots_L23_20180531(paths_EN,paths_NH,type,npSub)


       for K=1:length(paths_EN)
            [ROI_pos_plot_EN{K},responses_plot_EN{K}] = colorCodeByWhiskResponse_norm( paths_EN{K},type,npSub,[] );
            
        end
        
        ROI_pos_plot_EN=cat(1,ROI_pos_plot_EN{:});
        responses_plot_EN=cat(2,responses_plot_EN{:});
        
        % all ROIs
        for K=1:length(paths_NH)
            [ROI_pos_plot_NH{K},responses_plot_NH{K}] = colorCodeByWhiskResponse_norm( paths_NH{K},type,npSub,[] );
        end
        
        ROI_pos_plot_NH=cat(1,ROI_pos_plot_NH{:});
        responses_plot_NH=cat(2,responses_plot_NH{:});
        
 
%%

% 
% 
% figure; hold on
% subplot(1,2,1)
% colorCodeROIs_norm_mean(ROI_pos_plot_NH,'scale_var',responses_plot_NH,'type',type)
% title('NH')
% 
% subplot(1,2,2)
% colorCodeROIs_norm_mean(ROI_pos_plot_EN,'scale_var',responses_plot_EN,'type',type)
% title('EN')
% 
%% kmeans cluster points, voronoi diagram - response size vs. um distance

pool = gcp;                      % Invokes workers
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);

nClust_EN=floor(length(responses_plot_EN)/10);
nClust_NH=floor(length(responses_plot_NH)/10);
numBins_EN=15;
numBins_NH=15;
% nClust_EN=500;
% nClust_NH=500;
% numBins=15;

% ROI_pos_plot_NH=fliplr(ROI_pos_plot_NH);
% ROI_pos_plot_EN=fliplr(ROI_pos_plot_EN);

%find cmap bounds
[pos_idx,ROIpos_clust_NH]=kmeans(ROI_pos_plot_NH,nClust_NH,'Options',options,'MaxIter',10000,...
    'Display','final','Replicates',10);
for i=1:nClust_NH
    Rs=responses_plot_NH(pos_idx==i);
    cVar_clust_NH(i)=mean(Rs);
end

[pos_idx,ROIpos_clust_EN]=kmeans(ROI_pos_plot_EN,nClust_EN,'Options',options,'MaxIter',10000,...
    'Display','final','Replicates',10);
for i=1:nClust_EN
    Rs=responses_plot_EN(pos_idx==i);
    cVar_clust_EN(i)=mean(Rs);
end


cmapBounds(1)=min([cVar_clust_EN,cVar_clust_NH]);
cmapBounds(2)=max([cVar_clust_EN,cVar_clust_NH]);
cmap_EN=brewermap(numBins_EN,'YlGnBu');
cmap_NH=brewermap(numBins_NH,'YlGnBu');
figure; hold on
subplot(1,2,1)
voronoi_ROI_pos(ROIpos_clust_NH,cVar_clust_NH,'numBins',numBins_NH,'cmapBounds',cmapBounds,'cmap',cmap_NH);
tmp=viscircles([0 0],150,'EdgeColor','k','LineStyle','--','LineWidth',2);
tmp.Children(1).Marker='none';
title('NH')

subplot(1,2,2)
voronoi_ROI_pos(ROIpos_clust_EN,cVar_clust_EN,'numBins',numBins_EN,'cmapBounds',cmapBounds,'cmap',cmap_EN);
tmp2=viscircles([0 0],150,'EdgeColor','k','LineStyle','--','LineWidth',2);
tmp2.Children(1).Marker='none';
title('EN')

end