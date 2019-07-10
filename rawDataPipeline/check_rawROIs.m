function  [filtTimeSeries_new,npfiltTimeSeries_new,ROIs_exclude]=check_rawROIs( filtTimeSeries,npfiltTimeSeries )
%GUI that loops through all ROIs, plotting raw fluorecence time series for
%ROI mask, neuropil mask, and ROI fluorescence post neuropil subtraction.
%black is ROI, red is mask, blue is subtraction
[ npNormTimeSeries ] = npSubtract_preDF( filtTimeSeries,npfiltTimeSeries,0.3 );

fns=fieldnames(filtTimeSeries);
cellNames=fieldnames(filtTimeSeries.(fns{1}));
ROIs_all={};
NPs_all={};
subs_all={};

for K=1:length(cellNames)
    ROI_all=[];
    NP_all=[];
    sub_all=[];
    for i=1:length(fns)
        ROI_all=[ROI_all,filtTimeSeries.(fns{i}).(cellNames{K})(10:end-10)];
        NP_all=[NP_all,npfiltTimeSeries.(fns{i}).(cellNames{K})(10:end-10)];
        sub_all=[sub_all,npNormTimeSeries.(fns{i}).(cellNames{K})(10:end-10)];
        
    end
    ROIs_all{K}=ROI_all;
    NPs_all{K}=NP_all;
    subs_all{K}=sub_all;
    
end

j=1;
inds_ex=zeros(1,length(cellNames));
g=figure; hold on
plot(1:length(ROIs_all{j}),ROIs_all{j},'k')
plot(1:length(NPs_all{j}),NPs_all{j}+1,'r')
plot(1:length(subs_all{j}),subs_all{j}+2,'b')




% Create push button
keep_btn = uicontrol('Style', 'pushbutton', 'String', 'keep',...
    'Units','normalized','Position', [0.55 0.005 0.045 0.05],...
    'Callback', @keep_ROI);

ex_btn = uicontrol('Style', 'pushbutton', 'String', 'exclude',...
    'Units','normalized','Position', [0.45 0.005 0.045 0.05],...
    'Callback', @exclude_ROI);






% Pause until figure is closed ---------------------------------------%
waitfor(g);
inds_ex=logical(inds_ex);
ROIs_exclude=cellNames(inds_ex);
for k=1:length(fns)
    filtTimeSeries_new.(fns{k})=rmfield(filtTimeSeries.(fns{k}),ROIs_exclude);
    npfiltTimeSeries_new.(fns{k})=rmfield(npfiltTimeSeries.(fns{k}),ROIs_exclude);
end

    function keep_ROI(source,event)
        j=j+1;
        if j<=length(cellNames)
            
            figure(g); hold off;
            plot(1:length(ROIs_all{j}),ROIs_all{j},'k')
            hold on
            plot(1:length(NPs_all{j}),NPs_all{j}+1,'r')
            plot(1:length(subs_all{j}),subs_all{j}+2,'b')
        end
    end

    function exclude_ROI(source,event)
        if j<=length(cellNames)
            inds_ex(j)=1;
        end
   
        j=j+1;
        if j<=length(cellNames)
            
            figure(g); hold off;
            plot(1:length(ROIs_all{j}),ROIs_all{j},'k')
            hold on
            plot(1:length(NPs_all{j}),NPs_all{j}+1,'r')
            plot(1:length(subs_all{j}),subs_all{j}+2,'b')
        end
             
    end
end
