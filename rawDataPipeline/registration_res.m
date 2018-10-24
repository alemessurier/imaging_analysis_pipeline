function [ex_image]=registration_res(varargin)
%% registration_res runs dftregistration on sets of tiff stacks collected with ScanImage,
% saves registered stacks.
%
% INPUTs:
%       dir_to_reg -     directory to look for raw tiffs (required, must be
%                        listed before dir_to_save, set in analysisTemplate as 'dir_raw')
%
%       dir_to_save -    directory to save registered tiffs (required, must be
%                        listed after dir_to_reg in list of inputs; set in 
%                        analysisTemplate as 'dir_reg')
%       
%       ref_stack_path - file path of stack to use as reference (optional;
%                        if not set, the reference stack will be the middle
%                        stack found in dir_to_reg)
%
%
%% function body

indsDir=cellfun(@(x)isdir(x),varargin,'Uni',1);
dirsReg=varargin(indsDir);
dir_to_reg=dirsReg{1};
dir_to_save=dirsReg{2};

indsRef=cellfun(@(x)strfind(x,'.tif'),varargin,'Uni',0);
indsRef=cellfun(@(x)~isempty(x),indsRef,'Uni',1);

cd(dir_to_reg)
imFiles=dir('*.tif');
stackNames=arrayfun(@(x)x.name,imFiles,'Uni',0);

if sum(indsRef)>0
    ref_stack_path=varargin{indsRef};
else
    ref_stack_path=strcat(dir_to_reg,stackNames{ceil(length(stackNames)/2)});
end

ref_stack = LoadTIFF_SI5(ref_stack_path);
ref=mean(ref_stack,3);
ex_image=ref;
imagesc(ref);

fft2Ref=fft2(ref);
for K=1:length(stackNames);
    display(K)
    FileTif=strcat(dir_to_reg,stackNames{K});
    nameToWrite=strcat('dft',stackNames{K});
    [curr_stack,Metadata]=LoadTIFF_SI5(FileTif);
    [Metadata] = ReadTIFFHeader_SI5(Metadata);
    reg_stack=zeros(Metadata.RowPixels,Metadata.ColPixels,Metadata.NumFrames);
    reg_stack=uint16(reg_stack);
    cd(dir_to_save);
    error=zeros(1,size(curr_stack,3));
    parfor i=1:size(curr_stack,3)
        curr_im=curr_stack(:,:,i);
        [dft_out,reg]=dftregistration(fft2Ref,fft2(curr_im),100);
        reg_im=abs(ifft2(reg));
        reg_stack(:,:,i)=reg_im;
        error(i)=dft_out(1);
    end
    Metadata.dftRegError=error;
    WriteTIFF(reg_stack,Metadata,strcat(dir_to_save,nameToWrite))
end
end
