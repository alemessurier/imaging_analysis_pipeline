function [ neuropil_mask,ROI_positions ] = HW_makeNpMask( ROI_positions,ring,gap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% remove empty ROI masks
indsEmpty=sum(sum(ROI_positions));
indsEmpty=indsEmpty(:);
indsEmpty=indsEmpty==0;
ROI_positions=ROI_positions(:,:,~indsEmpty);

pixOfGap=gap;
%pixel of the width of npmask ring
pixOfRim=ring;
%pixel of the width of gap between ROI and npmask 
countMatrix=sum(ROI_positions,3);
allROIMatrix=logical(countMatrix);
%obtain the logic of all ROIs, to make sure npmask won't go into any ROIs

for i=1:size(ROI_positions,3)
findROIPos=regionprops(ROI_positions(:,:,i),'BoundingBox'); % rectagular box encasing the ellipse of ROI
ROIPos=findROIPos.BoundingBox;
maskMarginOut=[ROIPos(1)-(pixOfGap+pixOfRim) ROIPos(2)-(pixOfGap+pixOfRim)...
    ROIPos(3)+2*(pixOfGap+pixOfRim) ROIPos(4)+2*(pixOfGap+pixOfRim)];
maskMarginIn=[ROIPos(1)-pixOfGap ROIPos(2)-pixOfGap ...
    ROIPos(3)+2*pixOfGap ROIPos(4)+2*pixOfGap];
hTmp=imshow(ones(512));
%set(gcf,'visible','off')
h=get(hTmp,'Parent');
neuropil_Rect_Out=imellipse(h,maskMarginOut);
neuropil_Rect_In=imellipse(h,maskMarginIn);
neuropil_Vert_Out=getVertices(neuropil_Rect_Out);
neuropil_Vert_In=getVertices(neuropil_Rect_In);
neuropil_mask_Out=poly2mask(neuropil_Vert_Out(:,1),neuropil_Vert_Out(:,2),512,512);
neuropil_mask_In=poly2mask(neuropil_Vert_In(:,1),neuropil_Vert_In(:,2),512,512);
neuropil_masktmp=neuropil_mask_Out-neuropil_mask_In;
neuropil_mask(:,:,i)=neuropil_masktmp & ~allROIMatrix;
end

% remove overlapping neuropil mask
% numNP=size(neuropil_mask,3);
% countNPMatrix=sum(neuropil_mask,3);
% overlap = countNPMatrix > 1;
% neuropil_mask(repmat(overlap, 1, 1, numNP)) = 0;

end

