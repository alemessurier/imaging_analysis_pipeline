function [ rawTimeSeries] = get_fluoTimeSeries_OL( varargin )
% Extracts raw fluorescence time series from ROIs, and subtracts neuropil
% signal surrouding each ROI
% 
% INPUTS: imPath - path to raw movie
% 
%         ROI_positions - cell array of ROI positions (output of
%         'define_ROIsAML')
%
% OUTPUS: rawTimeSeries - data structure containing raw fluorescence time
%         series, organized by movie name and ROI name.
%         
% Amy LeMessurier 2015

if numel(varargin)==2;
    ROI_positions=varargin{1};
    rawStack=varargin{2};
    xpix=size(rawStack,2);
    ypix=size(rawStack,1);
else
    [fname,imPath] = uigetfile('.tif');
    imPath=strcat(imPath,fname);
    [rawStack, Metadata] = LoadTIFF_SI5(imPath);
    xpix=Metadata.RowPixels;
    ypix=Metadata.ColPixels;
    
  

    if numel(varargin)==0
        [fname,ROIPath] = uigetfile('.mat');
        load(strcat(ROIPath,fname),'ROI_positions');
    else
        ROI_positions=varargin{1};
    end
    
end

for r=1:size(ROI_positions,3)

    roi_mask=logical(ROI_positions(:,:,r));
    ROI=strcat('ROI',num2str(r));
    
    for f=1:size(rawStack,3)
        frame=rawStack(:,:,f);
        roiFluor(f)=mean(frame(roi_mask));
    end
    rawTimeSeries.(ROI)=roiFluor;
    clear roiFluor
    
end

%      function [ neuropil_mask ] = make_neuropilMask( ROIvertices )
%         %UNTITLED4 Summary of this function goes here
%         %   Detailed explanation goes here
%         ROIpos=findEllipsePos(ROIvertices);
%         maskPosOut=[ROIpos(1)-30 ROIpos(2)-30 ROIpos(3)+60 ROIpos(4)+60];
%         maskPosIn=[ROIpos(1)-10 ROIpos(2)-10 ROIpos(3)+20 ROIpos(4)+20];
% %          rectangle('Position',maskPos,'Curvature',[1 1],'edgeColor','b','LineWidth',2);
% %         figure
%          hTmp=imshow(ones(512));
%         set(gcf,'visible','off')
%         h=get(hTmp,'Parent');
%         neuropil_Rect_Out=imellipse(h,maskPosOut);
%         neuropil_Rect_In=imellipse(h,maskPosIn);
%         neuropil_V_Out=getVertices(neuropil_Rect_Out);
%         neuropil_V_In=getVertices(neuropil_Rect_In);
%         neuropil_mask_Out=poly2mask(neuropil_V_Out(:,1),neuropil_V_Out(:,2),512,512);
%         neuropil_mask_In=poly2mask(neuropil_V_In(:,1),neuropil_V_In(:,2),512,512);
% %         ROI_mask=poly2mask(ROIvertices(:,1),ROIvertices(:,2),512,512);
%         neuropil_mask=neuropil_mask_Out-neuropil_mask_In;
%         neuropil_mask=logical(neuropil_mask);
% %         avg_im(neuropil_mask)=1;
% %         figure
% %         imshow(avg_im)
%     end


    
  
end

