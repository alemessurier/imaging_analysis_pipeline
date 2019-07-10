function varargout=get_NPcorr_mask(pathName,corrThreshold)


% Load in ROI masks, NP masks, and raw time series for both
cd(pathName);
ROI_path=dir('ROI_positions*');

if length(ROI_path)==1
    load([pathName,ROI_path.name])
    load([pathName,'ROI_timeSeries.mat'])
    
    
    % Read path to processed data from analysis_template
    cd(pathName);
    name=dir('analysisTemplate*');
    filecont=fileread(strcat(pathName,name.name));
    expr = '[^\n]*dir_processed=[^\n]*';
    dp_string = regexp(filecont,expr,'match');
    eval(dp_string{:});
    
    % make sure all processed movie files exist
    cd(dir_processed)
    allTifs=dir('*.tif');
    
    % make sure rawTimeSeries has all ROIs from ROI_positions
    fns=fieldnames(rawTimeSeries);
    cellNames=fieldnames(rawTimeSeries.(fns{1}));
    
    try
        % if # fns doesn't match # tifs, or # cellNames doesn't match #ROI
        % or NP masks, end function, give error message
        assert(length(cellNames)==size(ROI_positions,3),...
            'maskLoading:wrongSizeROI_pos',...
            'the number of ROIs doesnt match number of ROI masks for %s',pathName)
        assert(length(cellNames)==size(npMasks_basic,3),...%changed from npMasks_basic 20180814
            'maskLoading:wrongSizeNPmasks',...
            'the number of ROIs doesnt match number of NP masks for %s',pathName)
        assert(length(fns)==length(allTifs),...
            'numTifs:wrongNumMovies',...
            '# movies in dir_processed ~= lenghth of fns in rawTimeSeries for %s',pathName)
        
        % if no assertions fail, carry out function body
        curr_im=LoadTIFF_SI5([dir_processed,allTifs(1).name]);
        np_ROIcorr=calc_npMaskCorr_byPix(curr_im,npMasks_basic,cellNames,corrThreshold,rawTimeSeries.(fns{1}));
        
        for fn=2:length(fns) % loop through all movies
            curr_im=LoadTIFF_SI5([dir_processed,allTifs(fn).name]);
            np_ROIcorr=calc_npMaskCorr_byPix(curr_im,np_ROIcorr,cellNames,corrThreshold,rawTimeSeries.(fns{fn}));
        end
        varargout{1}=np_ROIcorr;
    catch ME
        varargout{1}=ME;
    end
else
    ME = MException('maskLoading:wrongNumFiles', ...
        'incorrect number of ROI_positions files for %s',pathName);
    varargout{1}=ME;
end



end

function [neuropil_mask_ROIcorr]=calc_npMaskCorr_byPix(curr_im,np_ROIcorr,cellNames,corrThreshold,rawTimeSeries)

cn_inds=1:length(cellNames);
neuropil_mask_ROIcorr=zeros(512,512,length(cellNames));

parfor i=1:size(np_ROIcorr,3)
    thisNPMask=np_ROIcorr(:,:,i); % pick one npMask at a time
    [x,y]=find(thisNPMask); % find index for every pixels in a npMask
    
    % reorder cellNames to start with this ROI
    cn_order=[cellNames{i};cellNames(cn_inds(cn_inds~=i))];
    
    for j=1:sum(sum(thisNPMask)) % for every pixels in one npMask
        thisPixelFluor=double(curr_im(x(j),y(j),:));
        
        for c=1:length(cellNames) % comparing the trace of one pixel to all ROI's traces
            cn=cn_order{c};
            [R,~]=corrcoef(thisPixelFluor,rawTimeSeries.(cn)); % R and p value of correlation
            
            if R(1,2)>corrThreshold
                thisNPMask(x(j),y(j))=0; % remove that pixel from npMask
                break
            end
            
        end
        
    end
    neuropil_mask_ROIcorr(:,:,i)=thisNPMask;
end

end