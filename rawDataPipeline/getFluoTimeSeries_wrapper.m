function [ rawTimeSeries,Metadata ] = getFluoTimeSeries_wrapper( dir_processed,  ROI_positions )
%Extracts raw fluorescecence of each ROI on each frame, saves to structure
%'rawTimeSeries' with fields for each ROI and each movie. extracts
%metadata for each movie from tiff header. 
cd(dir_processed)
imFiles=dir('*.tif');
stackNames=arrayfun(@(x)x.name,imFiles,'Uni',0);


parfor K=1:length(stackNames) % change parfor=>for if not parallelizing
    thisStack=stackNames{K};
    fns{K}=thisStack(1:(strfind(thisStack,'.')-1));
    [curr_im,tmpMetadata]=LoadTIFF_SI5(strcat(dir_processed,thisStack));
     [Metadata_all(K)] = ReadTIFFHeader_SI5(tmpMetadata);
    [ rawTimeSeries_all(K)] = get_fluoTimeSeries_OL(ROI_positions,curr_im);

end

for J=1:length(fns)
%     fName=['m',fns{J}];
    Metadata.(fns{J}) =Metadata_all(J);
    rawTimeSeries.(fns{J}) = rawTimeSeries_all(J);
end
    
% save(strcat(dir_reduced,'ROI_timeSeries.mat'),'rawTimeSeries')
% save(strcat(dir_reduced,'Metadata.mat'),'Metadata')

end

