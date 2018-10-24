function [ROI_positions,neuropil_mask]=label_ROIs(dir_ROIs,dir_ROI_positions,ring,gap,skip)
% Loops through movies, calls define_ROIsAML, and outputs ROIs
% for entire imaging session as ellipse positions.
%
% INPUTs:   dir_ROIs - directory containing registered 16bit tif stacks from
%                      one imaging field
%
%           dir_ROI_positions - (optional) directory containing existing
%                      ROI masks from a previous imaging session to adjust to this
%                      imaging field
%
%           ring - width (# pixels) to use for donut shaped neuropil masks surrounding
%                  each soma ROI
%
%           gap -  width (# pixels) of gap to leave between soma ROI and NP
%                  mask
%
%           skip - loop through every Nth movie, with N defined by skip
%
% OUTPUTs: ROI_positions - array containing all ROI masks, defined by
%                          binary 512x512 masks. Each mask is one soma.
%
%          neuropil_mask - array containing all NP masks, surrounding each 
%                          soma mask in ROI_positions.
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

    
for K=1:skip:length(stackNames);  % loop every N movie files
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







end