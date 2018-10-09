function pvals= make_CorrAnalysisPlots_NPsubCompare( paths_NH,paths_EN)


% for K=1:length(paths_NH)
%     [traceByStim,~,framesEvoked,permTestResults] = load_NPsub_data_L23( paths_NH{K},1 );
%     ROIs_NH{K}=find_sigROIs(permTestResults,traceByStim);
% end
% 
% for K=1:length(paths_EN)
%     [traceByStim,~,framesEvoked,permTestResults] = load_NPsub_data_L23( paths_EN{K},1 );
%     ROIs_EN{K}=find_sigROIs(permTestResults,traceByStim);
% end

ROIs_EN=[];
ROIs_NH=[];
for r=1:3
    
    [fullfield_NH(r),fullfield_EN(r),inBarr_NH(r),inBarr_EN(r),xBarr_NH(r),xBarr_EN(r)]=make_CorrAnalysisPlots( paths_NH,paths_EN,1,r,ROIs_NH,ROIs_EN );
%     corrPlots_functions(fullfield_NH(r),fullfield_EN(r),inBarr_NH(r),inBarr_EN(r),xBarr_NH(r),xBarr_EN(r),[],[],'L23')
end

%% in barrel NC
figure; hold on
h_NH=noiseCorrByDist_rval(inBarr_NH,15);
legend('r=0','r=0.3','r=0.7');
h_EN=noiseCorrByDist_rval(inBarr_EN,15);

h_NH(2).Marker='*';
h_NH(3).Marker='o';
for i=1:3
    h_EN(i).Color=[1 0 0];
end
h_EN(2).Marker='*';
h_EN(3).Marker='o';

xlabel('distance')
ylabel('noise corr coeff')
title('in barrel pairs')

%% in barrel SC
figure; hold on
h_NH=sigCorrByDist_rval(inBarr_NH,15);
legend('r=0','r=0.3','r=0.7');
h_EN=sigCorrByDist_rval(inBarr_EN,15);

h_NH(2).Marker='*';
h_NH(3).Marker='o';
for i=1:3
    h_EN(i).Color=[1 0 0];
end
h_EN(2).Marker='*';
h_EN(3).Marker='o';

xlabel('distance')
ylabel('signal corr coeff')
title('in barrel pairs')

%% x barrel NC
figure; hold on
h_NH=noiseCorrByDist_rval(xBarr_NH,15);
legend('r=0','r=0.3','r=0.7');
h_EN=noiseCorrByDist_rval(xBarr_EN,15);

h_NH(2).Marker='*';
h_NH(3).Marker='o';
for i=1:3
    h_EN(i).Color=[1 0 0];
end
h_EN(2).Marker='*';
h_EN(3).Marker='o';

xlabel('distance')
ylabel('noise corr coeff')
title('cross barrel pairs')

%% in barrel SC
figure; hold on
h_NH=sigCorrByDist_rval(xBarr_NH,15);
legend('r=0','r=0.3','r=0.7');
h_EN=sigCorrByDist_rval(xBarr_EN,15);

h_NH(2).Marker='*';
h_NH(3).Marker='o';
for i=1:3
    h_EN(i).Color=[1 0 0];
end
h_EN(2).Marker='*';
h_EN(3).Marker='o';

xlabel('distance')
ylabel('signal corr coeff')
title('cross barrel pairs')

%% compare means

for r=1:3
compare_4means_bar(inBarr_NH(r).SC,inBarr_EN(r).SC,xBarr_NH(r).SC,xBarr_EN(r).SC)
tmp=gca;
tmp.XTick=1:4;
tmp.XTickLabel={'NH, in col','EN, in col','NH, x col','EN, x col'};
tmp.XTickLabelRotation=45;
tmp.YLabel.String=['signal corr coeff, r=',num2str(r)];

figure; hold on
plot_4cdfs(inBarr_NH(r).SC,inBarr_EN(r).SC,xBarr_NH(r).SC,xBarr_EN(r).SC)
title(['signal corr coeff, r=',num2str(r)])
legend('NH, in col','EN, in col','NH, x col','EN, x col');


test=[inBarr_NH(r).SC; inBarr_EN(r).SC; xBarr_NH(r).SC; xBarr_EN(r).SC];
in_NH=repmat(1,size(inBarr_NH(r).SC));
in_EN=repmat(2,size(inBarr_EN(r).SC));
x_NH=repmat(3,size(xBarr_NH(r).SC));
x_EN=repmat(4,size(xBarr_EN(r).SC));
groupvar=[in_NH;in_EN;x_NH;x_EN];
[p,tbl,stats]=anova1(test,groupvar);
c=multcompare(stats)
pvals(r).SC=c;

% test different in means w/permutation test
pvals(r).SC_perm.ENinVx=permutationTest(inBarr_EN(r).SC,xBarr_EN(r).SC,10000);
pvals(r).SC_perm.NHinVx=permutationTest(inBarr_NH(r).SC,xBarr_NH(r).SC,10000);
pvals(r).SC_perm.ENvNH_in=permutationTest(inBarr_EN(r).SC,inBarr_NH(r).SC,10000);
pvals(r).SC_perm.ENvNH_x=permutationTest(xBarr_EN(r).SC,xBarr_NH(r).SC,10000);
end

for r=1:3
compare_4means_bar(inBarr_NH(r).NC,inBarr_EN(r).NC,xBarr_NH(r).NC,xBarr_EN(r).NC)
tmp=gca;
tmp.XTick=1:4;
tmp.XTickLabel={'NH, in col','EN, in col','NH, x col','EN, x col'};
tmp.XTickLabelRotation=45;
tmp.YLabel.String=['noise corr coeff, r=',num2str(r)];

figure; hold on
plot_4cdfs(inBarr_NH(r).NC,inBarr_EN(r).NC,xBarr_NH(r).NC,xBarr_EN(r).NC)
title(['noise corr coeff, r=',num2str(r)])
legend('NH, in col','EN, in col','NH, x col','EN, x col');

test=[inBarr_NH(r).NC; inBarr_EN(r).NC; xBarr_NH(r).NC; xBarr_EN(r).NC];
in_NH=repmat(1,size(inBarr_NH(r).NC));
in_EN=repmat(2,size(inBarr_EN(r).NC));
x_NH=repmat(3,size(xBarr_NH(r).NC));
x_EN=repmat(4,size(xBarr_EN(r).NC));
groupvar=[in_NH;in_EN;x_NH;x_EN];
[p,tbl,stats]=anova1(test,groupvar);
c=multcompare(stats)
pvals(r).NC=c;


pvals(r).NC_perm.ENinVx=permutationTest(inBarr_EN(r).NC,xBarr_EN(r).NC,10000);
pvals(r).NC_perm.NHinVx=permutationTest(inBarr_NH(r).NC,xBarr_NH(r).NC,10000);
pvals(r).NC_perm.ENvNH_in=permutationTest(inBarr_EN(r).NC,inBarr_NH(r).NC,10000);
pvals(r).NC_perm.ENvNH_x=permutationTest(xBarr_EN(r).NC,xBarr_NH(r).NC,10000);
end

%%
    function h=noiseCorrByDist_rval(group,binEdges)
        dists=cell(numel(group));
        NCs=dists;
        numROIs=dists;
        means=dists;
        sems=dists;
        
        for r=1:numel(group)
            [ dists{r},means{r},sems{r},numROIs{r},binEdges,NCs{r} ] = binVarByDist( group(r).dist,group(r).NC,binEdges );
            
            h(r)= errorbar(dists{r},means{r},sems{r},'k.-','LineWidth',1);
        
        end
        
         
    end

function h=sigCorrByDist_rval(group,binEdges)
        dists=cell(numel(group));
        SCs=dists;
        numROIs=dists;
        means=dists;
        sems=dists;
        
        for r=1:numel(group)
            [ dists{r},means{r},sems{r},numROIs{r},binEdges,SCs{r} ] = binVarByDist( group(r).dist,group(r).SC,binEdges );
            
            h(r)= errorbar(dists{r},means{r},sems{r},'k.-','LineWidth',1);
        
        end
        
         
    end

end
