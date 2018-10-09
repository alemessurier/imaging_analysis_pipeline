function [Metadata] = ReadTIFFHeader_SI5(Metadata)

% ReadTIFFHeader_SI5() extracts SI5 imaging information from the TIFF fileheader
% and stores it in new fields within the Metadata structure.  These fields are:
% FrameNumberRaw(frame)      Native frame number during acquisition, ie before any online frame averaging
% FrameNumberTIFF(frame)     Frame number in TIFF, after frame averaging
% FrameTimeInMovie(frame)    Time (sec) since start-of-movie trigger.  Used for stimulus alignment.
% acqNumAveragedFrames       Number of frames averaged online during acquisition before saving TIFF
% acqScanFramePeriod         Interval (s) between raw single frame
% acqZoomFactor              Integer 1 to n.  Used to determine scale bar.
% MovieChannels              Vector containing each imaging channel interleaved frame-by-frame in the TIFF.
%                            Note: number of channels interleaved in TIFF file is length(MovieChannels).

% You must first have run LoadTIFF_SI5 to read the Metadata.

NumberImages = length(Metadata.InfoImage);
Metadata.FrameNumberRaw = zeros([NumberImages 1]);       % 1-d vectors
Metadata.FrameNumberTIFF = zeros([NumberImages 1]);
Metadata.FrameTimeInMovie = zeros([NumberImages 1]);
Metadata.MovieChannels = zeros(4);          % 
str = Metadata.InfoImage(1).ImageDescription;

% Movie-specific parameters, read once from the first header
expression = sprintf('acqNumAveragedFrames.+\n');
matchStr = regexp(str,expression,'match','dotexceptnewline');    % this returns whole line
tempindex = strfind(matchStr{1,1},'=');    % index of equal sign
Metadata.acqNumAveragedFrames = str2double(matchStr{1,1}(tempindex+1:end));     

expression = sprintf('scanFramePeriod.+\n');
matchStr = regexp(str,expression,'match','dotexceptnewline');    % this returns whole line
tempindex = strfind(matchStr{1,1},'=');    % index of equal sign
Metadata.acqScanFramePeriod = str2double(matchStr{1,1}(tempindex+1:end));     

expression = sprintf('zoomFactor.+\n');
matchStr = regexp(str,expression,'match','dotexceptnewline');    % this returns whole line
tempindex = strfind(matchStr{1,1},'=');    % index of equal sign
Metadata.acqZoomFactor = str2double(matchStr{1,1}(tempindex+1:end));     

expression = sprintf('channelsSave.+\n');
matchStr = regexp(str,expression,'match','dotexceptnewline');    % this returns whole line
tempindex = strfind(matchStr{1,1},'=');    % index of equal sign
matchStr = matchStr{1,1}(tempindex+2:end);   % contains [1] or [1;2] or [1;2;3;4] etc
Metadata.MovieChannels = eval(matchStr);   % this converts '[1;2]' string notation to column vector containing [1;2]
% Note on format of MovieChannels vector:
% The number of imaging channels in TIFF movie is length(Metadata.MovieChannels).  
% Imaging channels are stored in TIFF file in alternating frame, with frame 1 
% representating the imaging channel MovieChannels(1), frame 2 representing the 
% imaging channel MovieChannels(2), etc.

% Frame-specific parameters, read separately from header of each frame.
for frameindex = 1:NumberImages
    
    Metadata.FrameNumberTIFF(frameindex) = frameindex;      
    
    str = Metadata.InfoImage(frameindex).ImageDescription;
    
    expression = sprintf('Frame Timestamp.+\n');
    matchStr = regexp(str,expression,'match','dotexceptnewline');    % this returns whole line
    tempindex = strfind(matchStr{1,1},'=');    % index of equal sign
    Metadata.FrameTimeInMovie(frameindex) = str2double(matchStr{1,1}(tempindex+1:end));     
    
    expression = sprintf('Frame Number.+\n');
    matchStr = regexp(str,expression,'match','dotexceptnewline');    % this returns whole line
    tempindex = strfind(matchStr{1,1},'=');    % index of equal sign
    Metadata.FrameNumberRaw(frameindex) = str2double(matchStr{1,1}(tempindex+1:end));     
      
end

end

%% ScanImage5 TIFF Header Structure  

%    This is the structure of SI5-specific information stored in the TIFF Image header.    
%    All info is contained in InfoImage(TIFF_frame_number).ImageDescription
%    This field is a string with subfields separated by /n (=char(10))
%
% Frame Number =                5
% Frame Timestamp(s) =               0.133179812
% Acq Trigger Timestamp(s) =              -0.000083300
% Next File Marker Timestamp(s) =    230584300921.369320000
% DC Overvoltage = 0
% scanimage.SI5.VERSION_MAJOR = 5
% scanimage.SI5.VERSION_MINOR = 0
% scanimage.SI5.acqBeamOverScan = 0
% scanimage.SI5.acqMode = 'grab'
% scanimage.SI5.acqModeStartTime = [2015 6 3 17 2 50.833]
% scanimage.SI5.acqNumAveragedFrames = 5
% scanimage.SI5.acqNumFrames = 3000
% scanimage.SI5.acqStartTime = [2015 6 3 17 2 52.583]
% scanimage.SI5.acqState = 'grab'
% scanimage.SI5.acqsPerLoop = 100
% scanimage.SI5.beamDirectMode = false
% scanimage.SI5.beamFlybackBlanking = true
% scanimage.SI5.beamLengthConstants = Inf
% scanimage.SI5.beamNumBeams = 1
% scanimage.SI5.beamPowers = 63
% scanimage.SI5.beamPzAdjust = 0
% scanimage.SI5.bidirectionalAcq = true
% scanimage.SI5.bscope2FlipperMirrorPosition = 'pmt'
% scanimage.SI5.bscope2GalvoGalvoMirrorInPath = 0
% scanimage.SI5.bscope2GalvoResonantMirrorInPath = 1
% scanimage.SI5.bscope2PmtGains = [0 0 0 0]
% scanimage.SI5.bscope2PmtPowersOn = [false false false false]
% scanimage.SI5.bscope2PmtTripped = [0 0 0 0]
% scanimage.SI5.bscope2PmtValsSet = false
% scanimage.SI5.bscope2ScanAlign = 0
% scanimage.SI5.cfgFilename = 'D:\CONFIG\ScanImage\SI 5.0\install.cfg'
% scanimage.SI5.chan1LUT = [-48 688]
% scanimage.SI5.chan2LUT = [0 8191]
% scanimage.SI5.chan3LUT = [0 32767]
% scanimage.SI5.chan4LUT = [0 32767]
% scanimage.SI5.channelOffsets = [19 243 0 0]
% scanimage.SI5.channelsAutoReadOffsets = true
% scanimage.SI5.channelsDisplay = [1;2]
% scanimage.SI5.channelsInputRange = {[-1 1] [-1 1] [-1 1] [-1 1]}
% scanimage.SI5.channelsSave = [1;2]
% scanimage.SI5.channelsSubtractOffset = [true true true true]
% scanimage.SI5.errorCondition = false
% scanimage.SI5.errorConditionIdentifiers = {}
% scanimage.SI5.errorConditionMessages = {}
% scanimage.SI5.fastZAcquisitionDelay = 0
% scanimage.SI5.fastZActive = false
% scanimage.SI5.fastZAllowLiveBeamAdjust = false
% scanimage.SI5.fastZDiscardFlybackFrames = false
% scanimage.SI5.fastZEnable = false
% scanimage.SI5.fastZFillFraction = []
% scanimage.SI5.fastZFramePeriodAdjustment = -100
% scanimage.SI5.fastZImageType = 'XY-Z'
% scanimage.SI5.fastZNumDiscardFrames = 0
% scanimage.SI5.fastZNumVolumes = 1
% scanimage.SI5.fastZPeriod = []
% scanimage.SI5.fastZScanType = 'sawtooth'
% scanimage.SI5.fastZSettlingTime = 0
% scanimage.SI5.fillFraction = 0.9
% scanimage.SI5.fillFractionTime = 0.712867
% scanimage.SI5.flybackLinesPerFrame = 16
% scanimage.SI5.focusDuration = Inf
% scanimage.SI5.frameAcqFcnDecimationFactor = 1
% scanimage.SI5.linesPerFrame = 512
% scanimage.SI5.loggingFramesPerFile = Inf
% scanimage.SI5.loopAcqInterval = 10
% scanimage.SI5.motorSecondMotorZEnable = false
% scanimage.SI5.numFrames = 1
% scanimage.SI5.overvoltageStatus = false
% scanimage.SI5.pixelsPerLine = 512
% scanimage.SI5.scanAngleMultiplierSlow = 1
% scanimage.SI5.scanFramePeriod = 0.0334177
% scanimage.SI5.scanPixelTimeMaxMinRatio = 3
% scanimage.SI5.scanPixelTimeMean = 8.84277e-08
% scanimage.SI5.scanPixelTimeStats.pixelTimeRatio = 3
% scanimage.SI5.scanPixelTimeStats.meanPixelTime = 8.84277e-08
% scanimage.SI5.scanShiftSlow = 0
% scanimage.SI5.shutterHoldOpenOnStack = false
% scanimage.SI5.stackNumSlices = 1
% scanimage.SI5.stackReturnHome = true
% scanimage.SI5.stackStartCentered = false
% scanimage.SI5.stackUseStartPower = false
% scanimage.SI5.stackUserOverrideLz = false
% scanimage.SI5.stackZStepSize = 1
% scanimage.SI5.triggerExternalEdges = {'rising' 'rising' 'rising'}
% scanimage.SI5.triggerExternalTerminals = {'PFI1' '' ''}
% scanimage.SI5.triggerTypeExternal = 1
% scanimage.SI5.userFunctionsCfg = <nonscalar struct/object>
% scanimage.SI5.userFunctionsOverride = <nonscalar struct/object>
% scanimage.SI5.userFunctionsUsr = <nonscalar struct/object>
% scanimage.SI5.usrFilename = 'D:\CONFIG\ScanImage\SI 5.0\install.usr'
% scanimage.SI5.zoomFactor = 3     




