function make_labeledMovie( path_movie,Stimuli,stimFramesAll,whisk,sampRate )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% index stimtimes by whisker identity
stimLength=ceil(0.5*sampRate(1));

stims=[whisk, 'blank'];



indsBS=strfind(path_movie,'\');
indsBS=indsBS(end);
inds=strfind(path_movie,'.');
fn=path_movie((indsBS+1):(inds-1));
dir_processed=path_movie(1:indsBS);
TIFF_Filename=strcat(dir_processed,fn,'_labeled.tif');

[currIm,Metadata]=LoadTIFF_SI5(path_movie);

filtMovie=zeros(size(currIm));

for t=3:(size(currIm,3)-2)
    filtMovie(:,:,t)=median(currIm(:,:,(t-2):(t+2)),3);
end
filtMovie=uint16(filtMovie);

stimOrder=Stimuli.(fn).Label+1;
stimFrames=stimFramesAll.(fn);
stimLabels=stims(stimOrder);


position=[20 60; 20 60];


for i=1:length(stimFrames)
    frames=stimFrames(i):(stimFrames(i)+stimLength);
    for s=1:length(frames)
        im=insertText(filtMovie(:,:,frames(s)),position,stimLabels(i),'TextColor','white','FontSize',36,'BoxOpacity',0);
        filtMovie(:,:,frames(s))=rgb2gray(im);
    end
    
    
    
    
end


WriteTIFF( filtMovie, Metadata, TIFF_Filename )






end





