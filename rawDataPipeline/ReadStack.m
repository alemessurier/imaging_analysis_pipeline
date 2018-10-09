function [ FinalImage ] = ReadStack(stack_to_read,tif_type)
%ReadStack reads a tiff stack into matlab
% Inputs: 
    % stack_to_Read: the path of the tif to be opened
    % tif_type: specify 'uint8' or 'uint16' for 8 or 16 bit, respectively
    %
    InfoImage=imfinfo(stack_to_read);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    
    FinalImage=zeros(nImage,mImage,NumberImages,tif_type);
    for i=1:NumberImages
        FinalImage(:,:,i)=imread(stack_to_read,'Index',i);
    end

end

