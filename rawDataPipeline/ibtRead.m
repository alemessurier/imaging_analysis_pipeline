function varargout=ibtRead(fn, showStruct)
% WARNING due to a serious bug in igor file encoding that prevents large files otherwise, this function assumes that the
% sweep duration DOES NOT CHANGE DURING RECORDING (does not change in
% number of samples).
%
% input fn - .ibt BRECHT file filename i.e. 'LFP.ibt'
% input showStruct -optional: 1= print out structure in command window,otherwise
%       print nothing/normally
% output: if nargout=1 variable
%           i.e. output=ibtRead('LFP01.ibt')
%           output is one struct containing all info and data
%         elseif nargout=2 variables
%           i.e. [data, param]=ibtRead('LFP01.ibt')
%           first variable is a struct with raw sweep data (data.raw)
%           second variable is a struct with all remaining parameter
%           information collected by BRECHT igor program and stored in .ibt
%         Note: assigning 2 output variables is useful for .mat
%         files, for example, the 'param' variable described above can be
%         loaded separately from the data: load(filename,'param')
%         This is much faster than loading all the data every time!
%
% For BRECHT_BI V2.01 and higher

% initialize some vectors:
if nargin<2
        showStruct=0;
end
fmStartFreq=[];fmStopFreq=[];filterCutoff=[];

% Brecht IBT read:
fid=fopen(fn);
fseek(fid,0,'bof');
magicNum1=fread(fid,1,'int16'); %21 = BRECHT file  (distinct from 11 = ECCLES)
this_wheader_posn= fread(fid,1,'uint32');

brecht.date=sprintf('%d_%d_%d, %d:%d:%d GMT',local_igorTime2vec(fread(fid,1,'single')));
ylabel=fread(fid,20,'*char')';
xlabel=fread(fid,20,'*char')';
recFN=fread(fid,30,'*char')';  %Filename at time of recording (changed size to 30. 12/12/14 AML)
brecht.fn=fn;
brecht.recordedFN=recFN(1:strfind(recFN,'|')-1);
brecht.xlabel=xlabel(1:strfind(xlabel,'|')-1);
brecht.ylabel=ylabel(1:strfind(ylabel,'|')-1);

%Wave Header
% if magicNum1==11 %Seek back to header position
    fseek(fid,this_wheader_posn,'bof')
% end
magicNum=fread(fid,1,'int16'); %2 bytes (should be equal to 12 for brechtjr 2.2 files)
sweepNum=fread(fid,1,'int16'); %sweeps indexed from 0
% sweepNum=0;

while ~isempty(magicNum)
    
    ind=sweepNum+1; %sweeps indexed from 1    
    sweepNumbers(ind)=sweepNum;
    numSamples(ind)=fread(fid,1,'single'); %4 bytes (32-bit floating point)
    scaleFactor(ind)=fread(fid,1,'uint32');
    adjAmpGain(ind)= fread(fid,1,'single');
    kHz(ind)= fread(fid,1,'single'); %sampling rate in kHz
    recordingMode(ind)=fread(fid,1,'single');
    dx(ind)=fread(fid,1,'single');
    sweepTime(ind)=fread(fid,1,'single'); %(seconds)
    
    %For Each piezo, Collect the following current injection info:
    % Not currently saved into output sruct 
    for i=1:3
        flag(i)= fread(fid,1,'int16'); %flag= fread(fid,1,'single')
        val(i)= fread(fid,1,'single');
        stimNum(i)=fread(fid,1,'single');
        trialNum(i)=fread(fid,1,'single');
    end
    
    DCflag=fread(fid,1,'single');
    DCval=fread(fid,1,'single');
    
    % 	// whisker stimulus & trial information (new to BRECHT)
    stim_num(ind)=fread(fid,1,'single');
    trial_number(ind)=fread(fid,1,'single');
    recordingDepth(ind)=fread(fid,1,'single'); %Recording Depth in microns
                                               %Save x and y coords too?
    sweepPoint(ind)=fread(fid,1,'uint32'); %
    %A kludge to get around the igor pointer overflow that occurs when recording
    %over ~33792000 samples
%     if ind==3
%         fseekDist=diff(sweepPoint); %Measure the sweep duration in samples
%     end
%     
%     if ind>=3
%         sweepPointtemp=fread(fid,1,'uint32'); %// pointer to recorded sweep
%         sweepPoint(ind)=sweepPoint(ind-1)+fseekDist;  %Literally add fseek distance.
%     else
%         sweepPoint(ind)=fread(fid,1,'uint32'); %
%     end
%     
    nextHeaderTemp=fread(fid,1,'uint32'); %// pointer to next waveheader
    prevHeaderTemp=fread(fid,1,'uint32'); %// pointer to previous waveheader
    
    %Stimulus information
    whiskerblock_magicnumber=fread(fid,1,'int16'); %should = 24
    if magicNum1==21 || magicNum1==22 %BRECHT only (as of 8-22-2011)
        piezoLabel{ind}=fread(fid,2,'*char')'; %There is a bug in igor that assumes that each stimulus has a different whisker label.....
        piezoNumber(ind)=fread(fid,1,'single');
        stimAmp(ind)=fread(fid,1,'single'); %(microns)
        stimDur(ind)=fread(fid,1,'single'); %ms
        stimShape(ind)=fread(fid,1,'single'); %microns
        stimRiseFall(ind)=fread(fid,1,'single'); %ms
        prestimDelay(ind)=fread(fid,1,'single'); %ms
        impulseNum(ind)=fread(fid,1,'single'); %largely unused by my code
        isi(ind)=fread(fid,1,'single');  %largely unused by my code
        if magicNum1==22 %Recording manipulator coordinates
            manipCoords(1:5,ind)=fread(fid,5,'single');
        end
    end
    
    %Go to sweep and capture it:
    fseek(fid,sweepPoint(ind),'bof');
    % display(ind) -- might be having trouble at this point with
    % calibration .ibt files
    magicNum=fread(fid,1,'int16');
    
    t=fread(fid,numSamples(ind),'int16');
    sweepData(ind,1:numSamples(ind))=t;
    
    %This is actually the beginning of each sweep header:
    magicNum=fread(fid,1,'int16'); %2 bytes
    sweepNum=fread(fid,1,'int16'); %sweeps indexed from 0
end

brecht.sweeps=sweepData;
brecht.sweepNum=sweepNumbers;
brecht.sweepTimes=sweepTime;
brecht.recDepth=recordingDepth;

if magicNum1==21 || magicNum1==22  %BRECHT
    brecht.stim.piezoLabel=piezoLabel;
    brecht.stim.piezoNumber=piezoNumber;
    brecht.stim.shape=stimShape;
    brecht.stim.num=stim_num;
    brecht.stim.preDelay=prestimDelay;
    brecht.stim.amp=stimAmp;
    brecht.stim.dur=stimDur;
    brecht.stim.riseFall=stimRiseFall;
    brecht.stim.numImpulses=impulseNum;
    brecht.stim.isi=isi;
    if magicNum1==22
        brecht.manip.xpos=manipCoords(1,:);
        brecht.manip.ypos=manipCoords(2,:);
        brecht.manip.zpos=manipCoords(3,:);
        brecht.manip.y0pos=manipCoords(4,:);
        brecht.manip.z0pos=manipCoords(5,:);
    end
else % Piezo calibration parameters for manual input of stim (ECCLES)
    brecht.stim.piezoLabel='2';
    brecht.stim.piezoNumber=2;
    brecht.stim.shape=4;
    brecht.stim.num=0;
    brecht.stim.preDelay=20;
    brecht.stim.amp=5;
    brecht.stim.dur=1000;
    brecht.stim.riseFall=0;
    brecht.stim.numImpulses=0;
    brecht.stim.isi=0;
    brecht.recordedFN=brecht.fn(max(findstr('/',brecht.fn))+1:end-4); %Sometimes this doesn't get recorded properly
end

%Assumes that these settings are not changing from trial to trial
brecht.settings.samples=numSamples(1);
brecht.settings.sweepDur=numSamples(1)/kHz(1)/1000;
brecht.settings.scaleFactor=scaleFactor(1);
brecht.settings.adjAmpGain=adjAmpGain(1);
brecht.settings.fskHz=kHz(1); %samling rate in kHz
brecht.settings.recordingMode=recordingMode(1);
brecht.settings.dx=dx(1);

fclose(fid);

if nargout==2
    data.raw=single(brecht.sweeps);  %No need to save these as doubles... right?
    data.fn=brecht.fn;
    varargout{1}=data;
    param=rmfield(brecht,'sweeps');
    varargout{2}=param;
else
    varargout={brecht};
end

if showStruct
    dispStruct(brecht)
end

%% Relevent Excerpt from BRECHT igor code:
% commandstr = "Redimension/N="+num2str(no_samples)+" savesweep"
% 	Execute commandstr
% 	
% 	scale_factor = 30000 / ADrange[channel]				// calculate scale factor for conversion to 16-bits (max value = 30,000 out of 65,536 total & 32,768 before sign change occurs)	
% 	savesweep = wavename*scale_factor					// scale sweep into 16-bit field format (0->30,000 out of usable 32,768) 
% 	
% 	commandstr = ExptFileName[channel][path]+".ibt"
% 	
% 	Open/A /P=savepath /T="IGT-" refnum commandstr
% 
% 	if (sweep_number== 0)			// First sweep is sweep # 0
% 	
% 		absolute_time = datetime
% 		
% 		FSetPos refnum, 0				
% 		FBinWrite/F=2 refnum, fheader_magicnumber				// ECCLES fileheader magicnumber = 11
% 		FStatus refnum
% 		backfillposn[channel][path] = V_filePos							// save byte address to wave0ptr
% 		variable temp = 0
% 		FBinWrite/U/F=3 refnum, temp						// wave0ptr (F=3 bytes unsigned)   (NULL)
% 		FBinWrite/F=4 refnum, absolute_time				// time of first file write (absolute)
% 		
% 		ydataname = "mV or pA"							// can change sweep - to -sweep, and will be deduced from channel mode on each sweep
% 		xdataname = "msec"
% 		exptname = ExptFileName[channel][path]												
% 		ydataname += separator; xdataname += separator; exptname += separator
% 		do
% 			ydataname += " "
% 		while (strlen(ydataname)<20)						// pad with spaces to 20 chars //
% 		do
% 			xdataname += " "
% 		while (strlen(xdataname)<20)						// pad with spaces to 20 chars //
% 		do
% 			exptname += " "								// pad with spaces to 20 chars //		
% 		while (strlen(exptname)<20)
% 		
% 		FBinWrite refnum, ydataname					// ydataname[20]			
% 		FBinWrite refnum, xdataname					// xdataname[20]			// UNUSED
% 		FBinWrite refnum, exptname					// Exptname[20]
% 		prev_wheader_posn[channel][path] = 0						// there is no previous wheader
% 	endif
% 												
% 	// write the wave header //
% 	FStatus refnum
% 	this_wheader_posn = V_filePos					// save byte address of this wheader
% 		
% 	FSetPos refnum, backfillposn[channel][path]						// go back and fill in 
% 	FBinWrite/U/F=3 refnum, this_wheader_posn		// previous pointer to this waveheader
% 	FSetPos refnum, this_wheader_posn				// return to start of waveheader
% 	
% 	FBinWrite/F=2 refnum, wheader_magicnumber		// magicnumber = 2 (2 bytes)							// ECCLES wheader with AO information = 12
% 	FBinWrite/F=2 refnum, sweep_number				// sweep number (2 bytes)
% 	FBinWrite/F=4 refnum, no_samples					//was 2 bytes, changed to 4 bytes
% 	FBinWrite/U/F=3 refnum, scale_factor				// scale factor for decoding data
% 	FBinWrite/F=4 refnum, adj_amplifier_gain			// adj_amplifier gain  = front-panel gain * amplifier scaling.  Needed for decoding the wave data.
% 	FBinWrite/F=4 refnum, kHz						// sample rate in kHz
% 	FBinWrite/F=4 refnum, recordingmode				// RECORDING MODE  0 = OFF;  1 = CURRENT CLAMP;  2= VOLTAGE CLAMP
% 	FBinWrite/F=4 refnum, dx							// dx for calculating x-axis
% 	FBinWrite/F=4 refnum, time_of_sweep				// time of sweep (in number of seconds)
% 
% 	pulse = 0				
% 	Wave cmdflag = $("command_pulse_flag"+num2str(channel))
% 	Wave cmdval = $("command_pulse_value"+num2str(channel))
% 	Wave cmdstart = $("command_pulse_start"+num2str(channel))
% 	Wave cmddur = $("command_pulse_duration"+num2str(channel))
% 	do												// write CommandWaveOut information for this sweep.
% 		if (command_enabled[channel] * channel_mode[channel])				// if command was disabled, no command was output, so write 0 to the file.
% 			flag = cmdflag[pulse]; val = cmdval[pulse]; start = cmdstart[pulse]; dur = cmddur[pulse]	
% 		else
% 			flag = 0; val = 0; start = 0; dur = 0
% 		endif
% 		FBinWrite/F=2 refnum, flag
% 		FBinWrite/F=4 refnum, val
% 		FBinWrite/F=4 refnum, start
% 		FBinWrite/F=4 refnum, dur
% 		pulse += 1
% 	while (pulse < 3)
% 
% 	flag = cmdflag[3];  val = cmdval[3]					// flag and value for DC current injection -- treated specially.  No start/stop fields saved.
% 	FBinWrite/F=4 refnum, flag
% 	FBinWrite/F=4 refnum, val
% 		
% 	val = 0											// whisker stimulus & trial information (new to BRECHT)
% 	FBinWrite/F=4 refnum, stim_number				// stimulus number
% 	FBinWrite/F=4 refnum, trial_number					// trial number
% 	FBinWrite/F=4 refnum, recording_depth				// recording depth
% 	
% 	FStatus refnum
% 	tempptr = V_filePos + 12 + 64						// sweep ptr will be (12 bytes + 64 byte whisker stimulus block) ahead of this position.  each F/3 write takes 4 bytes! // 	
% 	FBinWrite/U/F=3 refnum, tempptr					// ptr to wavedata for this sweep (current position + 12 bytes) 
% 	FStatus refnum
% 	backfillposn[channel][path] = V_filePos				// backfillposn marks ptr to next waveheader
% 	FBinWrite/U/F=3 refnum, temp	
% 	variable t2=prev_wheader_posn[channel][path]				// ptr to next waveheader (NULL)
% 	FBinWrite/U/F=3 refnum,t2						// ptr to previous waveheader (NULL for sweep 1)
% 	prev_wheader_posn[channel][path] = this_wheader_posn		// save byte address of this header for next one
% 	
% 	// write the whisker stimulus block 					// UNIQUE TO BRECHT.  Only found in fheadermagicnumber = 21 BRECHT files
% 	FBinWrite/F=2 refnum, whiskerblock_magicnumber	// = 24 (Brecht)       This entire block is 64 bytes of data (36 bytes of data plus 28 bytes reserved)
% 	piezostr = PiezoLabel[stim_number]+"  "
% 	piezostr = piezostr[0,1]							// take first 2 characters.  This ensures that length is always 2 bytes, not shorter or longer.
% 	FBinWrite refnum, piezostr							// 2-character piezo label
% 	val = stim_piezo[stim_number]
% 	FBinWrite/F=4 refnum, val							// piezo number
% 	val = stim_amplitude[stim_number]
% 	FBinWrite/F=4 refnum, val							// stim amplitude (microns)
% 	val = stim_duration[stim_number]
% 	FBinWrite/F=4 refnum, val							// stim duration (ms)
% 	val = stim_shape[stim_number]
% 	FBinWrite/F=4 refnum, val							// stim shape (microns)
% 	val = stim_risefalltime[stim_number]
% 	FBinWrite/F=4 refnum, val							// stim rise fall time
% 	FBinWrite/F=4 refnum, WStim_Delay				// prestimulus delay (ms)
% 	FBinWrite/F=4 refnum, WStim_NumImpulses		// number of impulses
% 	FBinWrite/F=4 refnum, WStim_ISI					// interpulse interval
% 	commandstr = "                            "				// 28 bytes reserved for future use
% 	FBinWrite refnum, commandstr
% 	
% 	// write the wave data itself //
% 	FBinWrite/F=2 refnum, sweep_magicnumber			// ECCLES sweepdata magicnumber = 13 
% 	FBinWrite/F=2 refnum, savesweep					// the sweep itself (2 bytes per sample)
% 	
% 	// close the file
% 	Close refnum
% 	Killwaves savesweep
% 	
% 	//print "Wrote ch # ", channel, "path # ", path, " sweep", sweep_number, "to file ", commandstrefnum, time_of_sweep				// time of sweep (in number of seconds)

function DV=local_igorTime2vec(it)
% igor time -> date vector a la Matlab (MATLAB CENTRAL)
% The date/time is store as seconds since midnight, January 1, 1904.
Nsec_4year = 126230400;
idate_1988 = 21*uint32(Nsec_4year); % # seconds between 1-jan-1904 & 1-jan-1988
it = double(it-idate_1988); % smaller numbers-> use ordinary numbers ("doubles")
M4 = floor(it/Nsec_4year);
YR = double(1988+4*M4);
it = double(it-M4*Nsec_4year);
last_it = it;
while 1, % subtract years as long as they fit
    s = etime([YR+1 1 1 0 0 0], [YR 1 1 0 0 0]);
    if it<s, break; end % we went too far
    YR = YR+1;
    it = it-s;
end
MNTH = 1;
while 1, % subtract months as long as they fit
    s = etime([YR MNTH+1 1 0 0 0], [YR MNTH 1 0 0 0]);
    if it<s, break; end % we went too far
    MNTH = MNTH+1;
    it = it-s;
end
DAY = 1;
while 1, % subtract days as long as they fit
    s = etime([YR MNTH DAY+1 0 0 0], [YR MNTH DAY 0 0 0]);
    if it<s, break; end % we went too far
    DAY = DAY+1;
    it = it-s;
end
HR = floor(it/60^2); it = it-HR*60^2;
MIN = floor(it/60); SEC = it-MIN*60;
DV = [YR, MNTH, DAY, HR, MIN, SEC];
if ~isequal(last_it, etime(DV, [1988+4*M4 1 1 0 0 0])),
    error('date error');
end
DV = [MNTH, DAY, YR, HR, MIN, SEC];


function dispStruct(Xname,X)
%--- Modified from rec_structdisp(Xname,X) on FileExchange -BI

if nargin<2
    X=Xname;
    Xname='structure';
end

%-- PARAMETERS (Edit this) --%

ARRAYMAXROWS = 10;
ARRAYMAXCOLS = 10;
ARRAYMAXELEMS = 30;
CELLMAXROWS = 10;
CELLMAXCOLS = 10;
CELLMAXELEMS = 30;
CELLRECURSIVE = true;
% Xname = 'structure';
%----- PARAMETERS END -------%

disp(sprintf('%s :',Xname))
disp(X)
fprintf('\n')

if isstruct(X)
    F = fieldnames(X);
    nsub = length(F);
    Y = cell(1,nsub);
    subnames = cell(1,nsub);
    for i=1:nsub
        f = F{i};
        Y{i} = X.(f);
        subnames{i} = [Xname '.' f];
    end
elseif CELLRECURSIVE && iscell(X)
    nsub = numel(X);
    s = size(X);
    Y = X(:);
    subnames = cell(1,nsub);
    for i=1:nsub
        inds = s;
        globind = i-1;
        for k=1:length(s)
            inds(k) = 1+mod(globind,s(k));
            globind = floor(globind/s(k));
        end
        subnames{i} = [Xname '{' num2str(inds,'%i,') '}'];
    end
end

for i=1:nsub
    a = Y{i};
    if isstruct(a)
        if length(a)==1
            dispStruct(subnames{i},a)
        else
            for k=1:length(a)
                dispStruct([subnames{i} '(' num2str(k) ')'],a(k))
            end
        end
    elseif iscell(a)
        if size(a,1)<=CELLMAXROWS && size(a,2)<=CELLMAXCOLS && numel(a)<=CELLMAXELEMS
            dispStruct(subnames{i},a)
        end
    elseif size(a,1)>1 && size(a,1)<=ARRAYMAXROWS && size(a,2)<=ARRAYMAXCOLS && numel(a)<=ARRAYMAXELEMS
        disp([subnames{i} ':'])
        disp(a)
    end
end
