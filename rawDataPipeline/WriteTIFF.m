function WriteTIFF( MovieData, Metadata, TIFF_Filename )
%WriteTIFF will write an elementary TIFF movie file containing MovieData with
%the TIFF header information from Metadata.  TIFF_Filename should contain
%the full path.

% This is a direct interface to libtiff
t = Tiff(TIFF_Filename,'w');        % opens a new TIFF object for writing

MovieData = uint16(MovieData);      % convert to 16-bit if not already

% overall header
% For info on TIFF tags, see http://www.mathworks.com/help/matlab/ref/tiffclass.html
tagstruct.ImageLength     = size(MovieData,1);
tagstruct.ImageWidth      = size(MovieData,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 8;      % as per ScanImage5 TIFF format
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
tagstruct.ImageDescription = Metadata.InfoImage(1).ImageDescription;     % Long string including FrameTime and ScanImage settings
t.setTag(tagstruct)     % set tags
t.write(MovieData(:,:,1));   % writes first frame plus the metadata

% write each subsequent frame in a separate TIFF directory with its own tags
for (frame=2:size(MovieData,3))         % number of frames
    t.writeDirectory();         % set up next subdirectory (for next frame) and make it current
    tagstruct.ImageLength     = size(MovieData,1);
    tagstruct.ImageWidth      = size(MovieData,2);
    tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample   = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip    = 8;      % as used in ScanImage5 TIFF format
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software        = 'MATLAB';
    tagstruct.ImageDescription = Metadata.InfoImage(frame).ImageDescription;     
                    % Struct containing ScanImage specific header
                    % information for each frame
    t.setTag(tagstruct)     % write the tag
    
    t.write(MovieData(:,:,frame));   % writes this one frame plus the metadata
end
t.close();


end

