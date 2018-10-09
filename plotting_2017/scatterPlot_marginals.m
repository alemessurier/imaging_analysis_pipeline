function [figName,ax,binEdgesX,binEdgesY] = scatterPlot_marginals( X,Y,binEdgesX,binEdgesY )
%TOTESCOMBOPLOT Plots NH and EN binned variable by distance, with righthand axis
%for number of ROIs per bin, plus vertical histogram 
%   INPUTs:     dist_NH/dist_EN:    distance from barrel center (or
%                                   pairwise distance)
%           
%               var_NH/var_EN:      variable binned by distances in first
%                                   input
%       
%               sem_NH/sem_EN:      standard error for var_NH/EN
%           
%               num_ROIs_EN/NH:     number of ROIs per bin

%% 
figName=figure; 



% make left-handed plot for var/distance
ax(1)=axes('XAxisLocation','bottom',...
    'Color','none',...
     'Position',[0.1    0.1    0.6    0.6]);
 hold on
plot_scatterRLine( X,Y )



 linkaxes(ax,'x')
 
% make rotated histogram of NC as separate plot

ax(2)=axes('Position',[0.8    0.1    0.1    0.6],...
    'Units','normalized');
 binMin=ax(1).YLim(1);
 binMax=ax(1).YLim(2);
 if isempty(binEdgesY)
    binEdgesY=binMin:(binMax-binMin)/50:binMax;
 end
 
 histY=histc(Y,binEdgesY);
histY=histY/sum(histY);
 handlebar_Y=barh(binEdgesY,histY);
handlebar_Y.FaceColor='none';
handlebar_Y.EdgeColor='k';
handlebar_Y.LineWidth=1.5;
ax(2)=gca;
ax(2).YLim=ax(1).YLim;
% ax(2).YTick=[];
ax(2).XLim=[0 max(histY)];
% ax(2).XTick=[0 max(histY)];
% ax(2).XTickLabel=[0 max(histY)];
ax(2).YAxisLocation='left';
ax(2).Box='off';

% make rotated histogram of NC as separate plot

ax(3)=axes('Position',[0.1 0.8 0.6 0.1],...
    'Units','normalized');
 binMin=ax(1).XLim(1);
 binMax=ax(1).XLim(2);
 if isempty(binEdgesX)
    binEdgesX=binMin:(binMax-binMin)/50:binMax;
 end
 histX=histc(X,binEdgesX);
histX=histX/sum(histX);
 handlebar_Y=bar(binEdgesX,histX);
handlebar_Y.FaceColor='none';
handlebar_Y.EdgeColor='k';
handlebar_Y.LineWidth=1.5;
ax(3)=gca;
ax(3).XLim=ax(1).XLim;
% ax(3).XTick=[];
ax(3).YLim=[0 max(histX)];
% ax(3).XTick=[0 max(histY)];
% ax(3).XTickLabel=[0 max(histY)];
ax(3).Box='off';

linkaxes([ax(1),ax(2)],'y')
linkaxes([ax(1),ax(3)],'x')

end

