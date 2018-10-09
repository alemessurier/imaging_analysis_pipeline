function [ROI_positions,neuropil_mask]=label_ROIs_HW_AML(dir_ROIs,dir_ROI_positions,ring,gap)
% Loops through every fourth movie, calls define_ROIsAML, and outputs ROIs
% for entire imaging session as ellipse positions.
%
% INPUTs: dir_ROIs - directory containing registered 16bit tif stacks from
% one imaging field
%
% OUTPUTs: ROI_positions - array containing all ROIs, defined by positions
% vectors (4 vertices of ellipse). Each z is one ROI
%
% Amy LeMessurier 2014

imFiles=dir(dir_ROIs);   % list of movies
for i=1:length(imFiles)
    stackNames{i}=imFiles(i).name;   %stackNames=movie names
    dirNames(i)=imFiles(i).isdir;    %logic if is directory
end
stackNames=stackNames(~dirNames);    %stackNames not including folders

if ~isempty(dir_ROI_positions)  %has existing ROI_positions to compare
    load(strcat(dir_ROI_positions,'ROI_positions'),'ROI_positions');%(filename, variables)
    temp_im=double(ReadStack(strcat(dir_ROI_positions,'ex_image.tif'),'uint16')); % (file name to read, tif filetype)
    tmpIm2=double(ReadStack(strcat(dir_ROIs,stackNames{floor(length(stackNames)/2)}),'uint16'));  %load the middle movie
    toReg_im=mean(tmpIm2,3);  % example of current image section to compare with old section
    [shift,reg_im]=dftregistration(fft2(toReg_im),fft2(temp_im),100); % register example from old ROIs and example of current ROIs
    if shift(3)>0
        ROI_positions= ROI_positions(1:end-round(shift(3)),:,:);
            ROI_positions=cat(1,zeros(abs(round(shift(3))),size(ROI_positions,2),size(ROI_positions,3)),ROI_positions); 
    else
        ROI_positions= ROI_positions(abs(round(shift(3)))+1:end,:,:);
            ROI_positions=cat(1,ROI_positions,zeros(abs(round(shift(3))),size(ROI_positions,2),size(ROI_positions,3)));
    end
 
    
    if shift(4)<0
        ROI_positions=ROI_positions(:,abs(round(shift(4)))+1:end,:);
            ROI_positions=cat(2,ROI_positions,zeros(size(ROI_positions,1),abs(round(shift(4))),size(ROI_positions,3)));
    else
         ROI_positions=ROI_positions(:,1:end-abs(round(shift(4))),:);
            ROI_positions=cat(2,zeros(size(ROI_positions,1),abs(round(shift(4))),size(ROI_positions,3)),ROI_positions);
    end
else
    ROI_positions=zeros(512,512);
end

    
for K=1:5:length(stackNames);  % loop every 5 movie files
    imPath=strcat(dir_ROIs,stackNames{K});
    [im, ~] = LoadTIFF_SI5(imPath);
    ROIs_new=define_ROI_masks(im,ROI_positions); % ROI(:,:,1)=zeros(512,512) in function
    if isempty(dir_ROI_positions)&& K==1
        ROI_positions=ROIs_new(:,:,2:end);
    else
        ROI_positions=ROIs_new;
    end

end

for i=1:size(ROI_positions,3)
        emptyROIInd(i)=isempty(find(ROI_positions(:,:,i)));
end
ROI_positions(:,:,emptyROIInd)=[];





%neuropil mask new
%create a donut-shaped npMask beginning gap pixels from the outer most of
%ROIs and has the width of donut ring = ring pixels

[neuropil_mask]=HW_makeNpMask(ROI_positions,ring,gap); 




%remove pixel which have highly correlated activity with its neighbors,
%creating npmask for corr=1
% for i=1:size(neuropil_mask,3)
%     neuropil_mask_corr(:,:,i)=neuropil_mask(:,:,i) & ~correlation_map;
% end






%neuropil mask old
% this convolution method creats a really big npMask with no gaps between
% npMask and ROIs, not used in analysis now
% npMasksOld=makeNeuropilMasks(ROI_positions);

function NeuropilMasks=makeNeuropilMasks(ROIMasks)

numROIs=size(ROIMasks,3);
countMatrix = sum(ROIMasks, 3);                    % determine number of ROIs found in each pixel
overlap = countMatrix > 1;                         % define regions where ROIs overlap
ROIMasks(repmat(overlap, 1, 1, numROIs)) = 0;      % remove regions of overlap from ROI masks

g = exp(-(-10:10).^2/2/2^2);
maskb = conv2(g,g,double(logical(countMatrix)),'same')>.15;                        % dilation for border region around ROIs
[xi,yi] = meshgrid(1:size(ROIMasks,1),1:size(ROIMasks,2));
for rindex = 1:numROIs
    if sum(sum(ROIMasks(:,:,rindex)))==0
         NeuropilMasks(:,:,rindex) = zeros(size(ROIMasks,1),size(ROIMasks,2)); 
    else
        centroid=regionprops(ROIMasks(:,:,rindex),'centroid');
        centroids(rindex,:)=centroid.Centroid; %centroid is structure array, and centroid corordinate is under field "Centroid"
        for neuropilrad = 40:5:100
            M = (xi-centroids(rindex,1)).^2+(yi-centroids(rindex,2)).^2 < neuropilrad^2;    % mask of pixels within the radius
            NeuropilMasks(:,:,rindex) = M.*~maskb;                                          % remove ROIs and border regions
            if nnz(NeuropilMasks(:,:,rindex)) > 4000
                break
            end
        end
    end
end
end




end