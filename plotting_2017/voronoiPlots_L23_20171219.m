function voronoiPlots_L23_20171219(paths_EN,paths_NH,cells)


switch cells
    case 'all'
        % all ROIs
        for K=1:length(paths_EN)
            [ROI_pos_plot_EN{K},responses_plot_EN{K}] = colorCodeByWhiskResponse_norm( paths_EN{K} );
        end
        
        ROI_pos_plot_EN=cat(1,ROI_pos_plot_EN{:});
        responses_plot_EN=cat(2,responses_plot_EN{:});
        
        % all ROIs
        for K=1:length(paths_NH)
            [ROI_pos_plot_NH{K},responses_plot_NH{K}] = colorCodeByWhiskResponse_norm( paths_NH{K} );
        end
        
        ROI_pos_plot_NH=cat(1,ROI_pos_plot_NH{:});
        responses_plot_NH=cat(2,responses_plot_NH{:});
        
    case 'CW'
        % CWtuned only
        for K=1:length(paths_EN)
            [ROI_pos_plot_EN{K},responses_plot_EN{K}] = colorCodeByWhiskResponse_norm_SWorCW( paths_EN{K},'CW' );
        end
        
        ROI_pos_plot_EN=cat(1,ROI_pos_plot_EN{:});
        responses_plot_EN=cat(2,responses_plot_EN{:});
        
        
        % CWtuned only
        for K=1:length(paths_NH)
            [ROI_pos_plot_NH{K},responses_plot_NH{K}] = colorCodeByWhiskResponse_norm_SWorCW( paths_NH{K},'CW' );
        end
        
        ROI_pos_plot_NH=cat(1,ROI_pos_plot_NH{:});
        responses_plot_NH=cat(2,responses_plot_NH{:});
    case 'SW'
        % SWtuned only
        for K=1:length(paths_EN)
            [ROI_pos_plot_EN{K},responses_plot_EN{K}] = colorCodeByWhiskResponse_norm_SWorCW( paths_EN{K},'SW' );
        end
        
        ROI_pos_plot_EN=cat(1,ROI_pos_plot_EN{:});
        responses_plot_EN=cat(2,responses_plot_EN{:});
        
        
        % SWtuned only
        for K=1:length(paths_NH)
            [ROI_pos_plot_NH{K},responses_plot_NH{K}] = colorCodeByWhiskResponse_norm_SWorCW( paths_NH{K},'SW' );
        end
        
        ROI_pos_plot_NH=cat(1,ROI_pos_plot_NH{:});
        responses_plot_NH=cat(2,responses_plot_NH{:});
        
    case 'MW'
        % SWtuned only
        for K=1:length(paths_EN)
            [ROI_pos_plot_EN{K},responses_plot_EN{K}] = colorCodeByWhiskResponse_norm_MW( paths_EN{K},'MW' );
        end
        
        ROI_pos_plot_EN=cat(1,ROI_pos_plot_EN{:});
        responses_plot_EN=cat(2,responses_plot_EN{:});
        
        
        % SWtuned only
        for K=1:length(paths_NH)
            [ROI_pos_plot_NH{K},responses_plot_NH{K}] = colorCodeByWhiskResponse_norm_MW( paths_NH{K},'MW' );
        end
        
        ROI_pos_plot_NH=cat(1,ROI_pos_plot_NH{:});
        responses_plot_NH=cat(2,responses_plot_NH{:});
        
    case 'OW'
        % SWtuned only
        for K=1:length(paths_EN)
            [ROI_pos_plot_EN{K},responses_plot_EN{K}] = colorCodeByWhiskResponse_norm_MW( paths_EN{K},'OW' );
        end
        
        ROI_pos_plot_EN=cat(1,ROI_pos_plot_EN{:});
        responses_plot_EN=cat(2,responses_plot_EN{:});
        
        
        % SWtuned only
        for K=1:length(paths_NH)
            [ROI_pos_plot_NH{K},responses_plot_NH{K}] = colorCodeByWhiskResponse_norm_MW( paths_NH{K},'OW' );
        end
        
        ROI_pos_plot_NH=cat(1,ROI_pos_plot_NH{:});
        responses_plot_NH=cat(2,responses_plot_NH{:});
end

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
    Rs=responses_plot_NH(pos_idx==i);
    cVar_clust_NH(i)=mean(Rs);
end

[pos_idx,ROIpos_clust_EN]=kmeans(ROI_pos_plot_EN,nClust,'Options',options,'MaxIter',10000,...
    'Display','final','Replicates',10);
for i=1:nClust
    Rs=responses_plot_EN(pos_idx==i);
    cVar_clust_EN(i)=mean(Rs);
end


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

end