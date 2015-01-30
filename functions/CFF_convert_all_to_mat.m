function [ALLfileinfo] = CFF_convert_all_to_mat(ALLfilename, varargin)
% function [ALLfileinfo] = CFF_convert_all_to_mat(ALLfilename, varargin)
%
% DESCRIPTION
%
% converts Kongsberg EM series binary .all or .wcd data files (ALL) to a
% Matlab format (MAT), conserving all information from the original as it
% is. 
%
% INPUT VARIABLES
%
% - ALLfilename: .all or .wcd file to parse
% - varargin{1} (optional): output MAT file. If unspecified, the MAT file
% is saved in same folder as input ALL file and bears the same name except
% for its extension.
%
% OUTPUT VARIABLES
%
% - ALLfileinfo (optional): structure for description of the contents of
% the all file. Fields are:
%   * ALLfilename
%   * filesize (in bytes)
%   * datagsizeformat (endianness of the datagram size field)
%   * datagramsformat (endianness of the datagrams)
%   * datagramtypenumber (for each datagram, SIMRAD datagram type in decimal)
%   * datagramtype (for each datagram, SIMRAD datagram type description)
%   * parsedstatus (for each datagram, 1 if datagram has been parsed, 0 if not)
%   * datagramcounter (for each datagram, counter of this type of datagram in the file)
%   * datagramsize (for each datagram, datagram size in bytes)
%   * synccounter (for each datagram, number of bytes founds between this datagram and the previous one, indicative of sync error) 
%
% RESEARCH NOTES
%
% - PU Status output datagram structure seems different to the datagram manual
% description. Find the good description.#edit 21aug2013: updated to Rev Q. Need to be checked though.# 
%
% - output ALLfileinfo is growing inside loops and makes the function
% taking much longer. Run without if not required.
%
% - code currently lists the EM model numbers supported as a test for sync.
% Add your model number in the list if it's not currently there. It would
% be better to remove this test and try to sync on ETX and Checksum
% instead.
% 
% - to print out GPS datagrams (GGA), type: cell2mat(EM_Position.PositionInputDatagramAsReceived')
%
% - attitude datagrams contain several values of attitude. to pad cell values to allow plot, type:
% for ii = 1:length(EM_Attitude.TypeOfDatagram)
%     EM_Attitude.TimeInMillisecondsSinceRecordStart{ii} = [EM_Attitude.TimeInMillisecondsSinceRecordStart{ii};nan(max(EM_Attitude.NumberOfEntries)-length(EM_Attitude.TimeInMillisecondsSinceRecordStart{ii}),1)];
%     EM_Attitude.SensorStatus{ii} = [EM_Attitude.SensorStatus{ii};nan(max(EM_Attitude.NumberOfEntries)-length(EM_Attitude.SensorStatus{ii}),1)];
%     EM_Attitude.Roll{ii} = [EM_Attitude.Roll{ii};nan(max(EM_Attitude.NumberOfEntries)-length(EM_Attitude.Roll{ii}),1)];
%     EM_Attitude.Pitch{ii} = [EM_Attitude.Pitch{ii};nan(max(EM_Attitude.NumberOfEntries)-length(EM_Attitude.Pitch{ii}),1)];
%     EM_Attitude.Heave{ii} = [EM_Attitude.Heave{ii};nan(max(EM_Attitude.NumberOfEntries)-length(EM_Attitude.Heave{ii}),1)];
%     EM_Attitude.Heading{ii} = [EM_Attitude.Heading{ii};nan(max(EM_Attitude.NumberOfEntries)-length(EM_Attitude.Heading{ii}),1)];
% end
% % example: figure; grid on; plot(cell2mat(EM_Attitude.Roll))
%
% - to show soundspeed profile (if existing), type: figure;plot(cell2mat(EM_SoundSpeedProfile.Depth)./100, cell2mat(EM_SoundSpeedProfile.SoundSpeed)./10); grid on
%
% NEW FEATURES
%
% - v1.2:
%   - first added to SVN repository
%   - test existence of output directory and make it if not
% - v1.1:
%   - PU Status output datagram updated to Rev Q (but not tested)
%   - XYZ 88, Seabed Image 89, Raw Range and Angle 78, Water Column
%   Datagram & Network Attitude Velocity Datagram 110 now supported. 
% - v1.0:
%   - NA
% - v0.4.2:
%   - added synccounter
%   - changed ALLfileinfo field dimensions from 1 x nbping to nbping x 1
%   - very small change on test for input variables, for consistency with convxtf2all
% - v0.4:
%   - improved comments and general code
%   - test for byte ordering
%   - some EM_* structure names have been changed
%   - synchronization
%   - ETX check
%   - optional output information file
%   - separated start and stop installation parameters
% - v0.3.1:
%   - optional output MAT file name
%   - optional input machineformat
% - v0.2:
%   - optional soundspeed profiles supported
%
%%% 
% Alex Schimel, Deakin University
% Version 1.2 (28-04-2014)
%%%

%% supported systems:
emNumberList = [300; 2045; 3000; 3002; 3020]; %2045 is 2040c

%% if no input MATfilename, use same path and filename as ALL file
if nargin<2
    MATfilename = [ALLfilename(1:end-3) 'mat'];
else
    MATfilename = varargin{1};
end

%% check that output folder exists or create it
[p,n,e] = fileparts(MATfilename);
if ~exist(p,'dir')
    mkdir(p);
end
    
%% Checking byte ordering
% - Luciano's all files are in 'b'
% - Erik's all file is in 'l'
% - my converted files are in 'b'
% - DataDistrib files are in 'b' but datagram size in 'l'! We need to
% separate the byte ordering tests for these two types.

% opening file
[fid, message] = fopen(ALLfilename, 'r');

% number of bytes in file
temp = fread(fid,inf,'uint8');
filesize = length(temp);
clear temp
fseek(fid,0,-1);

% reading data from first datagram
while 1

    % read in little endian
    point = ftell(fid);
    nbDatagL = fread(fid,1,'uint32','l'); % number of bytes in datagram
        if isempty(nbDatagL), error('.all file parsing synchronization failed'); end %file finished
    fseek(fid,point,-1); % come back
    nbDatagB = fread(fid,1,'uint32','b'); % number of bytes in datagram
    stxDatag = fread(fid,1,'uint8'); % STX (always H02)
    typeDatag = fread(fid,1,'uint8'); % SIMRAD type of datagram
    emNumberL = fread(fid,1,'uint16','l'); % EM Model Number
    fseek(fid,-2,0); % come back
    emNumberB = fread(fid,1,'uint16','b'); % EM Model Number

    % trying to read ETX
    if fseek(fid,point+4+nbDatagL-3,-1) + 1 
        etxDatagL = fread(fid,1,'uint8'); % ETX (always H03)
    else
        etxDatagL = NaN;
    end
    if fseek(fid,point+4+nbDatagB-3,-1) + 1
        etxDatagB = fread(fid,1,'uint8'); % ETX (always H03)
    else
        etxDatagB = NaN;
    end

    % testing need for synchronization
    synchronized =    (sum(emNumberL==emNumberList) || sum(emNumberB==emNumberList)) ... 
                    & (etxDatagB==3 || etxDatagL==3) ...
                    & stxDatag==2;               
    if synchronized
        break
    else
        % trying to re-synchronize...
        fseek(fid,point+1,-1);
        continue
    end
end

% test for the byte ordering of the datagram size field
if etxDatagL == 3
    datagsizeformat = 'l';
elseif etxDatagB == 3
    datagsizeformat = 'b';
end

% test for byte ordering of datagrams
if sum(emNumberL==emNumberList)
    datagramsformat = 'l';
elseif sum(emNumberB==emNumberList)
    datagramsformat = 'b';
end

fclose(fid);

clear emNumberL emNumberB fid nbDatagL nbDatagB stxDatag typeDatag point etxDatagL etxDatagB synchronized

% create ouptut info file if required
if nargout
    ALLfileinfo.ALLfilename = ALLfilename;
    ALLfileinfo.filesize = filesize;
    ALLfileinfo.datagsizeformat = datagsizeformat;
    ALLfileinfo.datagramsformat = datagramsformat;
end


%% Reopening file with the good byte ordering

[fid, message] = fopen(ALLfilename, 'r',datagramsformat);

% intitializing datagram counter
kk = 0;

% initializing synchronization counter: the number of bytes that needed to
% be passed before this datagram appeared
syncn = 0;

%% Reading datagrams
while 1

    % new datagram begins, start reading
    point = ftell(fid);
    nbDatag = fread(fid,1,'uint32',datagsizeformat); % number of bytes in datagram
        if isempty(nbDatag), break; end % file finished, leave the loop  
    stxDatag = fread(fid,1,'uint8'); % STX (always H02)
    typeDatag = fread(fid,1,'uint8'); % SIMRAD type of datagram
    emNumber = fread(fid,1,'uint16'); % EM Model Number
    
    % test for synchronization
    % to pass, first data reading must show that:
    % - the number of bytes in following datagram doesn't overshoot file
    % size
    % - STX must be equal to 2.
    % - the EM model number must be in the list showed at beginning
    if nbDatag>filesize || stxDatag~=2 || ~sum(emNumber==emNumberList)
        fseek(fid,point+1,-1); % re-synchronizing 1 byte
        syncn = syncn+1; % update counter
        continue
    end
    
    switch typeDatag

%         case 49 % PU STATUS OUTPUT (31H)
%             
%             % datagrams counter
%             if ~exist('EM_PUStatus'), ii=1; else ii=size(EM_PUStatus.TypeOfDatagram,2)+1; end
% 
%             % parsing
%             EM_PUStatus.STX(ii)                                    = stxDatag;
%             EM_PUStatus.TypeOfDatagram(ii)                         = typeDatag;
%             EM_PUStatus.EMModelNumber(ii)                          = emNumber;
%             EM_PUStatus.Date(ii)                                   = fread(fid,1,'uint32');
%             EM_PUStatus.TimeSinceMidnightInMilliseconds(ii)        = fread(fid,1,'uint32');
%             EM_PUStatus.StatusDatagramCounter(ii)                  = fread(fid,1,'uint16');
%             EM_PUStatus.SystemSerialNumber(ii)                     = fread(fid,1,'uint16');                          
%             EM_PUStatus.PingRate(ii)                               = fread(fid,1,'uint16'); 
%             EM_PUStatus.PingCounterOfLatestPing(ii)                = fread(fid,1,'uint16');
%             EM_PUStatus.DistanceBetweenSwath(ii)                   = fread(fid,1,'uint8');
%             EM_PUStatus.SensorInputStatusUDPPort2(ii)              = fread(fid,1,'uint32');
%             EM_PUStatus.SensorInputStatusSerialPort1(ii)           = fread(fid,1,'uint32');
%             EM_PUStatus.SensorInputStatusSerialPort2(ii)           = fread(fid,1,'uint32');
%             EM_PUStatus.SensorInputStatusSerialPort3(ii)           = fread(fid,1,'uint32');
%             EM_PUStatus.SensorInputStatusSerialPort4(ii)           = fread(fid,1,'uint32');
%             EM_PUStatus.PPSStatus(ii)                              = fread(fid,1,'int8');
%             EM_PUStatus.PositionStatus(ii)                         = fread(fid,1,'int8');
%             EM_PUStatus.AttitudeStatus(ii)                         = fread(fid,1,'int8');
%             EM_PUStatus.ClockStatus(ii)                            = fread(fid,1,'int8');
%             EM_PUStatus.HeadingStatus (ii)                         = fread(fid,1,'int8');
%             EM_PUStatus.PUStatus(ii)                               = fread(fid,1,'uint8');
%             EM_PUStatus.LastReceivedHeading(ii)                    = fread(fid,1,'uint16');  
%             EM_PUStatus.LastReceivedRoll(ii)                       = fread(fid,1,'int16');  
%             EM_PUStatus.LastReceivedPitch(ii)                      = fread(fid,1,'int16');  
%             EM_PUStatus.LastReceivedHeave(ii)                      = fread(fid,1,'int16'); 
%             EM_PUStatus.SoundSpeedAtTransducer(ii)                 = fread(fid,1,'uint16'); 
%             EM_PUStatus.LastReceivedDepth(ii)                      = fread(fid,1,'uint32'); 
%             EM_PUStatus.AlongShipVelocity(ii)                      = fread(fid,1,'int16');
%             EM_PUStatus.AttitudeVelocitySensor(ii)                 = fread(fid,1,'uint8');
%             EM_PUStatus.MammalProtectionRamp(ii)                   = fread(fid,1,'uint8');
%             EM_PUStatus.BackscatterAtObliqueAngle(ii)              = fread(fid,1,'int8');
%             EM_PUStatus.BackscatterAtNormalIncidence(ii)           = fread(fid,1,'int8');
%             EM_PUStatus.FixedGain(ii)                              = fread(fid,1,'int8');
%             EM_PUStatus.DepthToNormalIncidence(ii)                 = fread(fid,1,'uint8');
%             EM_PUStatus.RangeToNormalIncidence(ii)                 = fread(fid,1,'uint16');
%             EM_PUStatus.PortCoverage(ii)                           = fread(fid,1,'uint8');
%             EM_PUStatus.StarboardCoverage(ii)                      = fread(fid,1,'uint8');
%             EM_PUStatus.SoundSpeedAtTransducerFoundFromProfile(ii) = fread(fid,1,'uint16');
%             EM_PUStatus.YawStabilization(ii)                       = fread(fid,1,'int16');
%             EM_PUStatus.PortCoverageOrAcrossShipVelocity(ii)       = fread(fid,1,'int16');
%             EM_PUStatus.StarboardCoverageOrDownwardVelocity(ii)    = fread(fid,1,'int16');
%             EM_PUStatus.EM2040CPUtemp(ii)                          = fread(fid,1,'int8'); 
%             EM_PUStatus.ETX(ii)                                    = fread(fid,1,'uint8');
%             EM_PUStatus.CheckSum(ii)                               = fread(fid,1,'uint16');
% 
%             % ETX check
%             if EM_PUStatus.ETX(ii)~=3
%                 error('wrong ETX value (EM_PUStatus)');
%             end
%                         
%             % file information
%             if nargout
%                 kk = kk+1;
%                 ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
%                 ALLfileinfo.datagramtype(kk,1) = {'PU STATUS OUTPUT (49, 31H)'};
%                 ALLfileinfo.parsedstatus(kk,1) = 1;
%                 ALLfileinfo.datagramcounter(kk,1) = ii;
%                 ALLfileinfo.datagramsize(kk,1) = nbDatag;
%                 ALLfileinfo.synccounter(kk,1) = syncn;
%             end
            
            
        case 65 % ATTITUDE (41H)

            % datagrams counter
            if ~exist('EM_Attitude'), ii=1; else ii=size(EM_Attitude.TypeOfDatagram,2)+1; end

            % parsing
            EM_Attitude.NumberOfBytesInDatagram(ii)                = nbDatag;
            EM_Attitude.STX(ii)                                    = stxDatag;
            EM_Attitude.TypeOfDatagram(ii)                         = typeDatag;
            EM_Attitude.EMModelNumber(ii)                          = emNumber;
            EM_Attitude.Date(ii)                                   = fread(fid,1,'uint32');
            EM_Attitude.TimeSinceMidnightInMilliseconds(ii)        = fread(fid,1,'uint32');
            EM_Attitude.AttitudeCounter(ii)                        = fread(fid,1,'uint16');
            EM_Attitude.SystemSerialNumber(ii)                     = fread(fid,1,'uint16');          
            EM_Attitude.NumberOfEntries(ii)                        = fread(fid,1,'uint16'); %N
            % repeat cycle: N entries of 12 bits
                temp = ftell(fid);
                N = EM_Attitude.NumberOfEntries(ii) ;
                EM_Attitude.TimeInMillisecondsSinceRecordStart{ii} = fread(fid,N,'uint16',12-2);
                fseek(fid,temp+2,'bof'); % to next data type
                EM_Attitude.SensorStatus{ii}                       = fread(fid,N,'uint16',12-2);
                fseek(fid,temp+4,'bof'); % to next data type
                EM_Attitude.Roll{ii}                               = fread(fid,N,'int16',12-2);
                fseek(fid,temp+6,'bof'); % to next data type
                EM_Attitude.Pitch{ii}                              = fread(fid,N,'int16',12-2);
                fseek(fid,temp+8,'bof'); % to next data type
                EM_Attitude.Heave{ii}                              = fread(fid,N,'int16',12-2);
                fseek(fid,temp+10,'bof'); % to next data type
                EM_Attitude.Heading{ii}                            = fread(fid,N,'uint16',12-2);
                fseek(fid,2-12,'cof'); % we need to come back after last jump
            EM_Attitude.SensorSystemDescriptor(ii)                 = fread(fid,1,'uint8');
            EM_Attitude.ETX(ii)                                    = fread(fid,1,'uint8');
            EM_Attitude.CheckSum(ii)                               = fread(fid,1,'uint16');

            % ETX check
            if EM_Attitude.ETX(ii)~=3
                error('wrong ETX value (EM_Attitude)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'ATTITUDE (65, 41H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end

        case 67 % CLOCK (43H)

            % datagrams counter
            if ~exist('EM_Clock'), ii=1; else ii=size(EM_Clock.TypeOfDatagram,2)+1; end

            % parsing
            EM_Clock.NumberOfBytesInDatagram(ii)                          = nbDatag;
            EM_Clock.STX(ii)                                              = stxDatag;
            EM_Clock.TypeOfDatagram(ii)                                   = typeDatag;
            EM_Clock.EMModelNumber(ii)                                    = emNumber;
            EM_Clock.Date(ii)                                             = fread(fid,1,'uint32');
            EM_Clock.TimeSinceMidnightInMilliseconds(ii)                  = fread(fid,1,'uint32');
            EM_Clock.ClockCounter(ii)                                     = fread(fid,1,'uint16');
            EM_Clock.SystemSerialNumber(ii)                               = fread(fid,1,'uint16');
            EM_Clock.DateFromExternalClock(ii)                            = fread(fid,1,'uint32');
            EM_Clock.TimeSinceMidnightInMillisecondsFromExternalClock(ii) = fread(fid,1,'uint32');
            EM_Clock.OnePPSUse(ii)                                        = fread(fid,1,'uint8');
            EM_Clock.ETX(ii)                                              = fread(fid,1,'uint8');
            EM_Clock.CheckSum(ii)                                         = fread(fid,1,'uint16');
            
            % ETX check
            if EM_Clock.ETX(ii)~=3
                error('wrong ETX value (EM_Clock)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'CLOCK (67, 43H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
            
        case 68 % DEPTH DATAGRAM (44H)

            % datagrams counter
            if ~exist('EM_Depth'), ii=1; else ii=size(EM_Depth.TypeOfDatagram,2)+1; end

            % parsing
            EM_Depth.NumberOfBytesInDatagram(ii)           = nbDatag;
            EM_Depth.STX(ii)                               = stxDatag;
            EM_Depth.TypeOfDatagram(ii)                    = typeDatag;
            EM_Depth.EMModelNumber(ii)                     = emNumber;
            EM_Depth.Date(ii)                              = fread(fid,1,'uint32');
            EM_Depth.TimeSinceMidnightInMilliseconds(ii)   = fread(fid,1,'uint32');
            EM_Depth.PingCounter(ii)                       = fread(fid,1,'uint16');
            EM_Depth.SystemSerialNumber(ii)                = fread(fid,1,'uint16'); 
            EM_Depth.HeadingOfVessel(ii)                   = fread(fid,1,'uint16'); 
            EM_Depth.SoundSpeedAtTransducer(ii)            = fread(fid,1,'uint16'); 
            EM_Depth.TransmitTransducerDepth(ii)           = fread(fid,1,'uint16');
            EM_Depth.MaximumNumberOfBeamsPossible(ii)      = fread(fid,1,'uint8');
            EM_Depth.NumberOfValidBeams(ii)                = fread(fid,1,'uint8'); %N
            EM_Depth.ZResolution(ii)                       = fread(fid,1,'uint8');
            EM_Depth.XAndYResolution(ii)                   = fread(fid,1,'uint8');
            EM_Depth.SamplingRate(ii)                      = fread(fid,1,'uint16'); % OR: EM_Depth.DepthDifferenceBetweenSonarHeadsInTheEM3000D(ii) = fread(fid,1,'int16');
            % repeat cycle: N entries of 16 bits
                temp = ftell(fid);
                N = EM_Depth.NumberOfValidBeams(ii);
                EM_Depth.DepthZ{ii}                        = fread(fid,N,'int16',16-2); % OR 'uint16' for EM120 and EM300
                fseek(fid,temp+2,'bof'); % to next data type
                EM_Depth.AcrosstrackDistanceY{ii}          = fread(fid,N,'int16',16-2);
                fseek(fid,temp+4,'bof'); % to next data type
                EM_Depth.AlongtrackDistanceX{ii}           = fread(fid,N,'int16',16-2);
                fseek(fid,temp+6,'bof'); % to next data type                    
                EM_Depth.BeamDepressionAngle{ii}           = fread(fid,N,'int16',16-2);
                fseek(fid,temp+8,'bof'); % to next data type                    
                EM_Depth.BeamAzimuthAngle{ii}              = fread(fid,N,'uint16',16-2);
                fseek(fid,temp+10,'bof'); % to next data type                    
                EM_Depth.Range{ii}                         = fread(fid,N,'uint16',16-2);
                fseek(fid,temp+12,'bof'); % to next data type                    
                EM_Depth.QualityFactor{ii}                 = fread(fid,N,'uint8',16-1);
                fseek(fid,temp+13,'bof'); % to next data type                         
                EM_Depth.LengthOfDetectionWindow{ii}       = fread(fid,N,'uint8',16-1);
                fseek(fid,temp+14,'bof'); % to next data type                         
                EM_Depth.ReflectivityBS{ii}                = fread(fid,N,'int8',16-1);
                fseek(fid,temp+15,'bof'); % to next data type                         
                EM_Depth.BeamNumber{ii}                    = fread(fid,N,'uint8',16-1);
                fseek(fid,1-16,'cof'); % we need to come back after last jump
            EM_Depth.TransducerDepthOffsetMultiplier(ii) = fread(fid,1,'int8');
            EM_Depth.ETX(ii)                             = fread(fid,1,'uint8');
            EM_Depth.CheckSum(ii)                        = fread(fid,1,'uint16');
            
            % ETX check
            if EM_Depth.ETX(ii)~=3,
                error('wrong ETX value (EM_Depth)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'DEPTH (68, 44H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
            
%         case 70 % RAW RANGE AND BEAM ANGLE (F) (46H)
%             
%             % to write...
                        
        case 71 % SURFACE SOUND SPEED (47H)

            % datagrams counter
            if ~exist('EM_SurfaceSoundSpeed'), ii=1; else ii=size(EM_SurfaceSoundSpeed.TypeOfDatagram,2)+1; end

            % parsing
            EM_SurfaceSoundSpeed.NumberOfBytesInDatagram(ii)           = nbDatag;
            EM_SurfaceSoundSpeed.STX(ii)                               = stxDatag;
            EM_SurfaceSoundSpeed.TypeOfDatagram(ii)                    = typeDatag;
            EM_SurfaceSoundSpeed.EMModelNumber(ii)                     = emNumber;
            EM_SurfaceSoundSpeed.Date(ii)                              = fread(fid,1,'uint32');
            EM_SurfaceSoundSpeed.TimeSinceMidnightInMilliseconds(ii)   = fread(fid,1,'uint32');
            EM_SurfaceSoundSpeed.SoundSpeedCounter(ii)                 = fread(fid,1,'uint16');
            EM_SurfaceSoundSpeed.SystemSerialNumber(ii)                = fread(fid,1,'uint16');            
            EM_SurfaceSoundSpeed.NumberOfEntries(ii)                   = fread(fid,1,'uint16'); %N     
            % repeat cycle: N entries of 4 bits
                temp = ftell(fid);
                N = EM_SurfaceSoundSpeed.NumberOfEntries(ii);
                EM_SurfaceSoundSpeed.TimeInSecondsSinceRecordStart{ii} = fread(fid,N,'uint16',4-2);
                fseek(fid,temp+2,'bof'); % to next data type
                EM_SurfaceSoundSpeed.SoundSpeed{ii}                    = fread(fid,N,'uint16',4-2);
                fseek(fid,2-4,'cof'); % we need to come back after last jump
            EM_SurfaceSoundSpeed.Spare(ii)                             = fread(fid,1,'uint8');
            EM_SurfaceSoundSpeed.ETX(ii)                               = fread(fid,1,'uint8');
            EM_SurfaceSoundSpeed.CheckSum(ii)                          = fread(fid,1,'uint16');

            % ETX check
            if EM_SurfaceSoundSpeed.ETX(ii)~=3
                error('wrong ETX value (EM_SurfaceSoundSpeed)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'SURFACE SOUND SPEED (71, 47H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
                        
%         case 72 % HEADING (48H)
%             
%             % to write...

        case 73 % INSTALLATION PARAMETERS - START (49H)

            % datagrams counter
            if ~exist('EM_InstallationStart'), ii=1; else ii=size(EM_InstallationStart.TypeOfDatagram,2)+1; end

            % parsing
            EM_InstallationStart.NumberOfBytesInDatagram(ii)         = nbDatag;
            EM_InstallationStart.STX(ii)                             = stxDatag;
            EM_InstallationStart.TypeOfDatagram(ii)                  = typeDatag;
            EM_InstallationStart.EMModelNumber(ii)                   = emNumber;
            EM_InstallationStart.Date(ii)                            = fread(fid,1,'uint32');
            EM_InstallationStart.TimeSinceMidnightInMilliseconds(ii) = fread(fid,1,'uint32');
            EM_InstallationStart.SurveyLineNumber(ii)                = fread(fid,1,'uint16');
            EM_InstallationStart.SystemSerialNumber(ii)              = fread(fid,1,'uint16');
            EM_InstallationStart.SerialNumberOfSecondSonarHead(ii)   = fread(fid,1,'uint16');

            % 18 bytes of binary data already recorded and 3 more to come = 21.
            % but nbDatag will always be even thanks to SpareByte. so
            % nbDatag is 22 if there is no ASCII data and more if there is
            % ASCII data. read the rest as ASCII (including SpareByte) with
            % 1 byte for 1 character. 
            EM_InstallationStart.ASCIIData{ii}                       = fscanf(fid, '%c', nbDatag-21);

            EM_InstallationStart.ETX(ii)                             = fread(fid,1,'uint8');
            EM_InstallationStart.CheckSum(ii)                        = fread(fid,1,'uint16');

            % ETX check
            if EM_InstallationStart.ETX(ii)~=3
                error('wrong ETX value (EM_InstallationStart)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'INSTALLATION PARAMETERS - START (73, 49H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
            
        case 78 % RAW RANGE AND ANGLE 78 (4EH)

            % datagrams counter
            if ~exist('EM_RawRangeAngle78'), ii=1; else ii=size(EM_RawRangeAngle78.TypeOfDatagram,2)+1; end

            % parsing
            EM_RawRangeAngle78.NumberOfBytesInDatagram(ii)           = nbDatag;
            EM_RawRangeAngle78.STX(ii)                               = stxDatag;
            EM_RawRangeAngle78.TypeOfDatagram(ii)                    = typeDatag;
            EM_RawRangeAngle78.EMModelNumber(ii)                     = emNumber;
            EM_RawRangeAngle78.Date(ii)                              = fread(fid,1,'uint32');
            EM_RawRangeAngle78.TimeSinceMidnightInMilliseconds(ii)   = fread(fid,1,'uint32');
            EM_RawRangeAngle78.PingCounter(ii)                       = fread(fid,1,'uint16');
            EM_RawRangeAngle78.SystemSerialNumber(ii)                = fread(fid,1,'uint16');
            EM_RawRangeAngle78.SoundSpeedAtTransducer(ii)            = fread(fid,1,'uint16'); 
            EM_RawRangeAngle78.NumberOfTransmitSectors(ii)           = fread(fid,1,'uint16'); %Ntx 
            EM_RawRangeAngle78.NumberOfReceiverBeamsInDatagram(ii)   = fread(fid,1,'uint16'); %Nrx
            EM_RawRangeAngle78.NumberOfValidDetections(ii)           = fread(fid,1,'uint16');
            EM_RawRangeAngle78.SamplingFrequencyInHz(ii)             = fread(fid,1,'float32');
            EM_RawRangeAngle78.Dscale(ii)                            = fread(fid,1,'uint32');
            % repeat cycle #1: Ntx entries of 24 bits
                temp = ftell(fid);
                C = 24;
                Ntx = EM_RawRangeAngle78.NumberOfTransmitSectors(ii);
                EM_RawRangeAngle78.TiltAngle{ii}                     = fread(fid,Ntx,'int16',C-2);
                fseek(fid,temp+2,'bof'); % to next data type       
                EM_RawRangeAngle78.FocusRange{ii}                    = fread(fid,Ntx,'uint16',C-2);
                fseek(fid,temp+4,'bof'); % to next data type
                EM_RawRangeAngle78.SignalLength{ii}                  = fread(fid,Ntx,'float32',C-4);
                fseek(fid,temp+8,'bof'); % to next data type
                EM_RawRangeAngle78.SectorTransmitDelay{ii}           = fread(fid,Ntx,'float32',C-4);
                fseek(fid,temp+12,'bof'); % to next data type
                EM_RawRangeAngle78.CentreFrequency{ii}               = fread(fid,Ntx,'float32',C-4);
                fseek(fid,temp+16,'bof'); % to next data type
                EM_RawRangeAngle78.MeanAbsorptionCoeff{ii}           = fread(fid,Ntx,'uint16',C-2);
                fseek(fid,temp+18,'bof'); % to next data type
                EM_RawRangeAngle78.SignalWaveformIdentifier{ii}      = fread(fid,Ntx,'uint8',C-1);
                fseek(fid,temp+19,'bof'); % to next data type
                EM_RawRangeAngle78.TransmitSectorNumberTxArrayIndex{ii} = fread(fid,Ntx,'uint8',C-1);
                fseek(fid,temp+20,'bof'); % to next data type
                EM_RawRangeAngle78.SignalBandwidth{ii}               = fread(fid,Ntx,'float32',C-4);
                fseek(fid,4-C,'cof'); % we need to come back after last jump
            % repeat cycle #2: Nrx entries of 16 bits    
                temp = ftell(fid);
                C = 16;
                Nrx = EM_RawRangeAngle78.NumberOfReceiverBeamsInDatagram(ii);
                EM_RawRangeAngle78.BeamPointingAngle{ii}             = fread(fid,Nrx,'int16',C-2);
                fseek(fid,temp+2,'bof'); % to next data type         
                EM_RawRangeAngle78.TransmitSectorNumber{ii}          = fread(fid,Nrx,'uint8',C-1);
                fseek(fid,temp+3,'bof'); % to next data type    
                EM_RawRangeAngle78.DetectionInfo{ii}                 = fread(fid,Nrx,'uint8',C-1);
                fseek(fid,temp+4,'bof'); % to next data type    
                EM_RawRangeAngle78.DetectionWindowLength{ii}         = fread(fid,Nrx,'uint16',C-2);
                fseek(fid,temp+6,'bof'); % to next data type    
                EM_RawRangeAngle78.QualityFactor{ii}                 = fread(fid,Nrx,'uint8',C-1);
                fseek(fid,temp+7,'bof'); % to next data type    
                EM_RawRangeAngle78.Dcorr{ii}                         = fread(fid,Nrx,'int8',C-1);
                fseek(fid,temp+8,'bof'); % to next data type    
                EM_RawRangeAngle78.TwoWayTravelTime{ii}              = fread(fid,Nrx,'float32',C-4);
                fseek(fid,temp+12,'bof'); % to next data type    
                EM_RawRangeAngle78.ReflectivityBS{ii}                = fread(fid,Nrx,'int16',C-2);
                fseek(fid,temp+14,'bof'); % to next data type    
                EM_RawRangeAngle78.RealTimeCleaningInfo{ii}          = fread(fid,Nrx,'int8',C-1);
                fseek(fid,temp+15,'bof'); % to next data type    
                EM_RawRangeAngle78.Spare{ii}                         = fread(fid,Nrx,'uint8',C-1);
                fseek(fid,1-C,'cof'); % we need to come back after last jump
            EM_RawRangeAngle78.Spare2(ii)                            = fread(fid,1,'uint8');
            EM_RawRangeAngle78.ETX(ii)                               = fread(fid,1,'uint8');
            EM_RawRangeAngle78.CheckSum(ii)                          = fread(fid,1,'uint16');
            
            % ETX check
            if EM_RawRangeAngle78.ETX(ii)~=3,
                error('wrong ETX value (EM_RawRangeAngle78)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'RAW RANGE AND ANGLE 78 (78, 4EH)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end

%         case 79 % QUALITY FACTOR DATAGRAM 79 (4FH)
%             
%             % to write...

        case 80 % POSITION (50H)

            % datagrams counter
            if ~exist('EM_Position'), ii=1; else ii=size(EM_Position.TypeOfDatagram,2)+1; end

            % parsing
            EM_Position.NumberOfBytesInDatagram(ii)         = nbDatag;
            EM_Position.STX(ii)                             = stxDatag;
            EM_Position.TypeOfDatagram(ii)                  = typeDatag;
            EM_Position.EMModelNumber(ii)                   = emNumber;
            EM_Position.Date(ii)                            = fread(fid,1,'uint32');
            EM_Position.TimeSinceMidnightInMilliseconds(ii) = fread(fid,1,'uint32');
            EM_Position.PositionCounter(ii)                 = fread(fid,1,'uint16');
            EM_Position.SystemSerialNumber(ii)              = fread(fid,1,'uint16');                
            EM_Position.Latitude(ii)                        = fread(fid,1,'int32');
            EM_Position.Longitude(ii)                       = fread(fid,1,'int32');
            EM_Position.MeasureOfPositionFixQuality(ii)     = fread(fid,1,'uint16');
            EM_Position.SpeedOfVesselOverGround(ii)         = fread(fid,1,'uint16');
            EM_Position.CourseOfVesselOverGround(ii)        = fread(fid,1,'uint16');
            EM_Position.HeadingOfVessel(ii)                 = fread(fid,1,'uint16');
            EM_Position.PositionSystemDescriptor(ii)        = fread(fid,1,'uint8');
            EM_Position.NumberOfBytesInInputDatagram(ii)    = fread(fid,1,'uint8');

            % next data size is variable. 34 bits of binary data already
            % recorded and 3 more to come = 37. read the rest as ASCII
            % (including SpareByte)
            EM_Position.PositionInputDatagramAsReceived{ii} = fscanf(fid, '%c', nbDatag-37);

            EM_Position.ETX(ii)                             = fread(fid,1,'uint8');
            EM_Position.CheckSum(ii)                        = fread(fid,1,'uint16');

            % ETX check
            if EM_Position.ETX(ii)~=3
                error('wrong ETX value (EM_Position)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'POSITION (80, 50H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
            
            
        case 82 % RUNTIME PARAMETERS (52H)

            % datagrams counter
            if ~exist('EM_Runtime'), ii=1; else ii=size(EM_Runtime.TypeOfDatagram,2)+1; end

            % parsing
            EM_Runtime.NumberOfBytesInDatagram(ii)                 = nbDatag;
            EM_Runtime.STX(ii)                                     = stxDatag;
            EM_Runtime.TypeOfDatagram(ii)                          = typeDatag;
            EM_Runtime.EMModelNumber(ii)                           = emNumber;
            EM_Runtime.Date(ii)                                    = fread(fid,1,'uint32');
            EM_Runtime.TimeSinceMidnightInMilliseconds(ii)         = fread(fid,1,'uint32');
            EM_Runtime.PingCounter(ii)                             = fread(fid,1,'uint16');
            EM_Runtime.SystemSerialNumber(ii)                      = fread(fid,1,'uint16');
            EM_Runtime.OperatorStationStatus(ii)                   = fread(fid,1,'uint8');
            EM_Runtime.ProcessingUnitStatus(ii)                    = fread(fid,1,'uint8');
            EM_Runtime.BSPStatus(ii)                               = fread(fid,1,'uint8');
            EM_Runtime.SonarHeadStatus(ii)                         = fread(fid,1,'uint8');
            EM_Runtime.Mode(ii)                                    = fread(fid,1,'uint8');
            EM_Runtime.FilterIdentifier(ii)                        = fread(fid,1,'uint8');
            EM_Runtime.MinimumDepth(ii)                            = fread(fid,1,'uint16');
            EM_Runtime.MaximumDepth(ii)                            = fread(fid,1,'uint16');
            EM_Runtime.AbsorptionCoefficient(ii)                   = fread(fid,1,'uint16');
            EM_Runtime.TransmitPulseLength(ii)                     = fread(fid,1,'uint16');
            EM_Runtime.TransmitBeamwidth(ii)                       = fread(fid,1,'uint16');
            EM_Runtime.TransmitPowerReMaximum(ii)                  = fread(fid,1,'int8');
            EM_Runtime.ReceiveBeamwidth(ii)                        = fread(fid,1,'uint8');
            EM_Runtime.ReceiveBandwidth(ii)                        = fread(fid,1,'uint8');
            EM_Runtime.ReceiverFixedGainSetting(ii)                = fread(fid,1,'uint8'); % OR mode 2
            EM_Runtime.TVGLawCrossoverAngle(ii)                    = fread(fid,1,'uint8');
            EM_Runtime.SourceOfSoundSpeedAtTransducer(ii)          = fread(fid,1,'uint8');
            EM_Runtime.MaximumPortSwathWidth(ii)                   = fread(fid,1,'uint16');
            EM_Runtime.BeamSpacing(ii)                             = fread(fid,1,'uint8');
            EM_Runtime.MaximumPortCoverage(ii)                     = fread(fid,1,'uint8');
            EM_Runtime.YawAndPitchStabilizationMode(ii)            = fread(fid,1,'uint8');
            EM_Runtime.MaximumStarboardCoverage(ii)                = fread(fid,1,'uint8');
            EM_Runtime.MaximumStarboardSwathWidth(ii)              = fread(fid,1,'uint16');
            EM_Runtime.DurotongSpeed(ii)                           = fread(fid,1,'uint16'); % OR: EM_Runtime.TransmitAlongTilt(ii) = fread(fid,1,'int16');
            EM_Runtime.HiLoFrequencyAbsorptionCoefficientRatio(ii) = fread(fid,1,'uint8'); % OR filter identifier 2
            EM_Runtime.ETX(ii)                                     = fread(fid,1,'uint8');
            EM_Runtime.CheckSum(ii)                                = fread(fid,1,'uint16');

            % ETX check
            if EM_Runtime.ETX(ii)~=3,
                error('wrong ETX value (EM_Runtime)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'RUNTIME PARAMETERS (82, 52H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end

        case 83 % SEABED IMAGE DATAGRAM (53H)

            % datagrams counter
            if ~exist('EM_SeabedImage'), ii=1; else ii=size(EM_SeabedImage.TypeOfDatagram,2)+1; end

            % parsing
            EM_SeabedImage.NumberOfBytesInDatagram(ii)         = nbDatag;
            EM_SeabedImage.STX(ii)                             = stxDatag;
            EM_SeabedImage.TypeOfDatagram(ii)                  = typeDatag;
            EM_SeabedImage.EMModelNumber(ii)                   = emNumber;
            EM_SeabedImage.Date(ii)                            = fread(fid,1,'uint32');
            EM_SeabedImage.TimeSinceMidnightInMilliseconds(ii) = fread(fid,1,'uint32');
            EM_SeabedImage.PingCounter(ii)                     = fread(fid,1,'uint16');
            EM_SeabedImage.SystemSerialNumber(ii)              = fread(fid,1,'uint16');
            EM_SeabedImage.MeanAbsorptionCoefficient(ii)       = fread(fid,1,'uint16'); % 'this field had earlier definition'
            EM_SeabedImage.PulseLength(ii)                     = fread(fid,1,'uint16'); % 'this field had earlier definition'
            EM_SeabedImage.RangeToNormalIncidence(ii)          = fread(fid,1,'uint16'); 
            EM_SeabedImage.StartRangeSampleOfTVGRamp(ii)       = fread(fid,1,'uint16');
            EM_SeabedImage.StopRangeSampleOfTVGRamp(ii)        = fread(fid,1,'uint16');
            EM_SeabedImage.NormalIncidenceBS(ii)               = fread(fid,1,'int8'); %BSN
            EM_SeabedImage.ObliqueBS(ii)                       = fread(fid,1,'int8'); %BSO
            EM_SeabedImage.TxBeamwidth(ii)                     = fread(fid,1,'uint16');
            EM_SeabedImage.TVGLawCrossoverAngle(ii)            = fread(fid,1,'uint8');
            EM_SeabedImage.NumberOfValidBeams(ii)              = fread(fid,1,'uint8'); %N    
            % repeat cycle: N entries of 6 bits
                temp = ftell(fid);
                N = EM_SeabedImage.NumberOfValidBeams(ii);
                EM_SeabedImage.BeamIndexNumber{ii}             = fread(fid,N,'uint8',6-1);
                fseek(fid,temp+1,'bof'); % to next data type
                EM_SeabedImage.SortingDirection{ii}            = fread(fid,N,'int8',6-1);
                fseek(fid,temp+2,'bof'); % to next data type
                EM_SeabedImage.NumberOfSamplesPerBeam{ii}      = fread(fid,N,'uint16',6-2); %Ns
                fseek(fid,temp+4,'bof'); % to next data type                    
                EM_SeabedImage.CentreSampleNumber{ii}          = fread(fid,N,'uint16',6-2);
                fseek(fid,2-6,'cof'); % we need to come back after last jump  
            Ns = [EM_SeabedImage.NumberOfSamplesPerBeam{ii}];
            for jj = 1:length(Ns)
                EM_SeabedImage.SampleAmplitudes(ii).beam{jj}   = fread(fid,Ns(jj),'int8');
            end
            % "spare byte if required to get even length (always 0 if used)"
            if floor(sum(Ns)/2) == sum(Ns)/2
                % even so far, since ETX is 1 byte, add a spare here
                EM_SeabedImage.Data.SpareByte(ii)              = fread(fid,1,'uint8');
            else
                % odd so far, since ETX is 1 bytes, no spare
                EM_SeabedImage.Data.SpareByte(ii) = NaN;
            end
            EM_SeabedImage.ETX(ii)                             = fread(fid,1,'uint8');
            EM_SeabedImage.CheckSum(ii)                        = fread(fid,1,'uint16');

            % ETX check
            if EM_SeabedImage.ETX(ii)~=3
                error('wrong ETX value (EM_SeabedImage)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'SEABED IMAGE DATAGRAM (83, 53H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
            
        case 85 % SOUND SPEED PROFILE (55H)

            % datagrams counter
            if ~exist('EM_SoundSpeedProfile'), ii=1; else ii=size(EM_SoundSpeedProfile.TypeOfDatagram,2)+1; end

            % parsing
            EM_SoundSpeedProfile.NumberOfBytesInDatagram(ii)                           = nbDatag;
            EM_SoundSpeedProfile.STX(ii)                                               = stxDatag;
            EM_SoundSpeedProfile.TypeOfDatagram(ii)                                    = typeDatag;
            EM_SoundSpeedProfile.EMModelNumber(ii)                                     = emNumber;
            EM_SoundSpeedProfile.Date(ii)                                              = fread(fid,1,'uint32');
            EM_SoundSpeedProfile.TimeSinceMidnightInMilliseconds(ii)                   = fread(fid,1,'uint32');            
            EM_SoundSpeedProfile.ProfileCounter(ii)                                    = fread(fid,1,'uint16');
            EM_SoundSpeedProfile.SystemSerialNumber(ii)                                = fread(fid,1,'uint16');
            EM_SoundSpeedProfile.DateWhenProfileWasMade(ii)                            = fread(fid,1,'uint32');
            EM_SoundSpeedProfile.TimeSinceMidnightInMillisecondsWhenProfileWasMade(ii) = fread(fid,1,'uint32');
            EM_SoundSpeedProfile.NumberOfEntries(ii)                                   = fread(fid,1,'uint16'); %N
            EM_SoundSpeedProfile.DepthResolution(ii)                                   = fread(fid,1,'uint16');
            % repeat cycle: N entries of 8 bits
                temp = ftell(fid);
                N = EM_SoundSpeedProfile.NumberOfEntries(ii);
                EM_SoundSpeedProfile.Depth{ii}                                         = fread(fid,N,'uint32',8-4);
                fseek(fid,temp+4,'bof'); % to next data type
                EM_SoundSpeedProfile.SoundSpeed{ii}                                    = fread(fid,N,'uint32',8-4);
                fseek(fid,4-8,'cof'); % we need to come back after last jump
            EM_SoundSpeedProfile.SpareByte(ii)                                         = fread(fid,1,'uint8');
            EM_SoundSpeedProfile.ETX(ii)                                               = fread(fid,1,'uint8');
            EM_SoundSpeedProfile.CheckSum(ii)                                          = fread(fid,1,'uint16');           

            % ETX check
            if EM_SoundSpeedProfile.ETX(ii)~=3
                error('wrong ETX value (EM_SoundSpeedProfile)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'SOUND SPEED PROFILE (85, 55H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
                        
        case 88 % XYZ 88 (58H)
            
            % datagrams counter
            if ~exist('EM_XYZ88'), ii=1; else ii=size(EM_XYZ88.TypeOfDatagram,2)+1; end

            % parsing
            EM_XYZ88.NumberOfBytesInDatagram(ii)           = nbDatag;
            EM_XYZ88.STX(ii)                               = stxDatag;
            EM_XYZ88.TypeOfDatagram(ii)                    = typeDatag;
            EM_XYZ88.EMModelNumber(ii)                     = emNumber;
            EM_XYZ88.Date(ii)                              = fread(fid,1,'uint32');
            EM_XYZ88.TimeSinceMidnightInMilliseconds(ii)   = fread(fid,1,'uint32');
            EM_XYZ88.PingCounter(ii)                       = fread(fid,1,'uint16');
            EM_XYZ88.SystemSerialNumber(ii)                = fread(fid,1,'uint16'); 
            EM_XYZ88.HeadingOfVessel(ii)                   = fread(fid,1,'uint16'); 
            EM_XYZ88.SoundSpeedAtTransducer(ii)            = fread(fid,1,'uint16'); 
            EM_XYZ88.TransmitTransducerDepth(ii)           = fread(fid,1,'float32');
            EM_XYZ88.NumberOfBeamsInDatagram(ii)           = fread(fid,1,'uint16');
            EM_XYZ88.NumberOfValidDetections(ii)           = fread(fid,1,'uint16');
            EM_XYZ88.SamplingFrequencyInHz(ii)             = fread(fid,1,'float32');
            EM_XYZ88.ScanningInfo(ii)                      = fread(fid,1,'uint8');
            EM_XYZ88.Spare1(ii)                            = fread(fid,1,'uint8');
            EM_XYZ88.Spare2(ii)                            = fread(fid,1,'uint8');
            EM_XYZ88.Spare3(ii)                            = fread(fid,1,'uint8');
            % repeat cycle: N entries of 20 bits
                temp = ftell(fid);
                C = 20;
                N = EM_XYZ88.NumberOfBeamsInDatagram(ii);
                EM_XYZ88.DepthZ{ii}                        = fread(fid,N,'float32',C-4);
                fseek(fid,temp+4,'bof'); % to next data type
                EM_XYZ88.AcrosstrackDistanceY{ii}          = fread(fid,N,'float32',C-4);
                fseek(fid,temp+8,'bof'); % to next data type
                EM_XYZ88.AlongtrackDistanceX{ii}           = fread(fid,N,'float32',C-4);
                fseek(fid,temp+12,'bof'); % to next data type                    
                EM_XYZ88.DetectionWindowLength{ii}         = fread(fid,N,'uint16',C-2);
                fseek(fid,temp+14,'bof'); % to next data type                    
                EM_XYZ88.QualityFactor{ii}                 = fread(fid,N,'uint8',C-1);
                fseek(fid,temp+15,'bof'); % to next data type                    
                EM_XYZ88.BeamIncidenceAngleAdjustment{ii}  = fread(fid,N,'int8',C-1);
                fseek(fid,temp+16,'bof'); % to next data type                    
                EM_XYZ88.DetectionInformation{ii}          = fread(fid,N,'uint8',C-1);
                fseek(fid,temp+17,'bof'); % to next data type                         
                EM_XYZ88.RealTimeCleaningInformation{ii}   = fread(fid,N,'int8',C-1);
                fseek(fid,temp+18,'bof'); % to next data type                         
                EM_XYZ88.ReflectivityBS{ii}                = fread(fid,N,'int16',C-2);
                fseek(fid,2-C,'cof'); % we need to come back after last jump
            EM_XYZ88.Spare4(ii)                            = fread(fid,1,'uint8');
            EM_XYZ88.ETX(ii)                               = fread(fid,1,'uint8');
            EM_XYZ88.CheckSum(ii)                          = fread(fid,1,'uint16');
            
            % ETX check
            if EM_XYZ88.ETX(ii)~=3,
                error('wrong ETX value (EM_XYZ88)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'XYZ 88 (88, 58H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
          
        case 89 % SEABED IMAGE DATA 89 (59H)
     
            % datagrams counter
            if ~exist('EM_SeabedImage89'), ii=1; else ii=size(EM_SeabedImage89.TypeOfDatagram,2)+1; end

            % parsing
            EM_SeabedImage89.NumberOfBytesInDatagram(ii)         = nbDatag;
            EM_SeabedImage89.STX(ii)                             = stxDatag;
            EM_SeabedImage89.TypeOfDatagram(ii)                  = typeDatag;
            EM_SeabedImage89.EMModelNumber(ii)                   = emNumber;
            EM_SeabedImage89.Date(ii)                            = fread(fid,1,'uint32');
            EM_SeabedImage89.TimeSinceMidnightInMilliseconds(ii) = fread(fid,1,'uint32');
            EM_SeabedImage89.PingCounter(ii)                     = fread(fid,1,'uint16');
            EM_SeabedImage89.SystemSerialNumber(ii)              = fread(fid,1,'uint16');
            EM_SeabedImage89.SamplingFrequencyInHz(ii)           = fread(fid,1,'float32');
            EM_SeabedImage89.RangeToNormalIncidence(ii)          = fread(fid,1,'uint16'); 
            EM_SeabedImage89.NormalIncidenceBS(ii)               = fread(fid,1,'int16'); %BSN
            EM_SeabedImage89.ObliqueBS(ii)                       = fread(fid,1,'int16'); %BSO
            EM_SeabedImage89.TxBeamwidthAlong(ii)                = fread(fid,1,'uint16');
            EM_SeabedImage89.TVGLawCrossoverAngle(ii)            = fread(fid,1,'uint16');
            EM_SeabedImage89.NumberOfValidBeams(ii)              = fread(fid,1,'uint16');      
            % repeat cycle: N entries of 6 bits
                temp = ftell(fid);
                C = 6;
                N = EM_SeabedImage89.NumberOfValidBeams(ii);
                EM_SeabedImage89.SortingDirection{ii}            = fread(fid,N,'int8',C-1);
                fseek(fid,temp+1,'bof'); % to next data type
                EM_SeabedImage89.DetectionInfo{ii}               = fread(fid,N,'uint8',C-1);
                fseek(fid,temp+2,'bof'); % to next data type
                EM_SeabedImage89.NumberOfSamplesPerBeam{ii}      = fread(fid,N,'uint16',C-2); %Ns
                fseek(fid,temp+4,'bof'); % to next data type                    
                EM_SeabedImage89.CentreSampleNumber{ii}          = fread(fid,N,'uint16',C-2);
                fseek(fid,2-C,'cof'); % we need to come back after last jump  
            Ns = [EM_SeabedImage89.NumberOfSamplesPerBeam{ii}];
            for jj = 1:length(Ns)
                EM_SeabedImage89.SampleAmplitudes(ii).beam{jj}   = fread(fid,Ns(jj),'int16');
            end
            EM_SeabedImage89.Spare(ii)                           = fread(fid,1,'uint8');
            EM_SeabedImage89.ETX(ii)                             = fread(fid,1,'uint8');
            EM_SeabedImage89.CheckSum(ii)                        = fread(fid,1,'uint16');

            % ETX check
            if EM_SeabedImage89.ETX(ii)~=3
                error('wrong ETX value (EM_SeabedImage89)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'SEABED IMAGE DATA 89 (89, 59H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end

%         case 102 % RAW RANGE AND BEAM ANGLE (f) (66H)
%             
%             % to write...

        case 104 % DEPTH (PRESSURE) OR HEIGHT DATAGRAM (68H) 

            % datagrams counter
            if ~exist('EM_Height'), ii=1; else ii=size(EM_Height.TypeOfDatagram,2)+1; end

            % parsing
            EM_Height.NumberOfBytesInDatagram(ii)         = nbDatag;
            EM_Height.STX(ii)                             = stxDatag;
            EM_Height.TypeOfDatagram(ii)                  = typeDatag;
            EM_Height.EMModelNumber(ii)                   = emNumber;
            EM_Height.Date(ii)                            = fread(fid,1,'uint32');
            EM_Height.TimeSinceMidnightInMilliseconds(ii) = fread(fid,1,'uint32');            
            EM_Height.HeightCounter(ii)                   = fread(fid,1,'uint16');
            EM_Height.SystemSerialNumber(ii)              = fread(fid,1,'uint16');
            EM_Height.Height(ii)                          = fread(fid,1,'int32');
            EM_Height.HeigthType(ii)                      = fread(fid,1,'uint8');
            EM_Height.ETX(ii)                             = fread(fid,1,'uint8');
            EM_Height.CheckSum(ii)                        = fread(fid,1,'uint16');           
            
            % ETX check
            if EM_Height.ETX(ii)~=3
                error('wrong ETX value (EM_Height)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'DEPTH (PRESSURE) OR HEIGHT (104, 68H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end

            
        case 105 % INSTALLATION PARAMETERS -  STOP (69H)

            % datagrams counter
            if ~exist('EM_InstallationStop'), ii=1; else ii=size(EM_InstallationStop.TypeOfDatagram,2)+1; end

            % parsing
            EM_InstallationStop.NumberOfBytesInDatagram(ii)         = nbDatag;
            EM_InstallationStop.STX(ii)                             = stxDatag;
            EM_InstallationStop.TypeOfDatagram(ii)                  = typeDatag;
            EM_InstallationStop.EMModelNumber(ii)                   = emNumber;
            EM_InstallationStop.Date(ii)                            = fread(fid,1,'uint32');
            EM_InstallationStop.TimeSinceMidnightInMilliseconds(ii) = fread(fid,1,'uint32');
            EM_InstallationStop.SurveyLineNumber(ii)                = fread(fid,1,'uint16');
            EM_InstallationStop.SystemSerialNumber(ii)              = fread(fid,1,'uint16');
            EM_InstallationStop.SerialNumberOfSecondSonarHead(ii)   = fread(fid,1,'uint16');

            % 18 bytes of binary data already recorded and 3 more to come = 21.
            % but nbDatag will always be even thanks to SpareByte. so
            % nbDatag is 22 if there is no ASCII data and more if there is
            % ASCII data. read the rest as ASCII (including SpareByte) with
            % 1 byte for 1 character.            
            EM_InstallationStop.ASCIIData{ii}                       = fscanf(fid, '%c', nbDatag-21);

            EM_InstallationStop.ETX(ii)                             = fread(fid,1,'uint8');
            EM_InstallationStop.CheckSum(ii)                        = fread(fid,1,'uint16');

            % ETX check
            if EM_InstallationStop.ETX(ii)~=3
                error('wrong ETX value (EM_InstallationStop)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'INSTALLATION PARAMETERS -  STOP (105, 69H)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
            
        case 107 % WATER COLUMN DATAGRAM (6BH)
            
            % datagrams counter
            if ~exist('EM_WaterColumn'), ii=1; else ii=size(EM_WaterColumn.TypeOfDatagram,2)+1; end

            % parsing
            EM_WaterColumn.NumberOfBytesInDatagram(ii)           = nbDatag;
            EM_WaterColumn.STX(ii)                               = stxDatag;
            EM_WaterColumn.TypeOfDatagram(ii)                    = typeDatag;
            EM_WaterColumn.EMModelNumber(ii)                     = emNumber;
            EM_WaterColumn.Date(ii)                              = fread(fid,1,'uint32');
            EM_WaterColumn.TimeSinceMidnightInMilliseconds(ii)   = fread(fid,1,'uint32');
            EM_WaterColumn.PingCounter(ii)                       = fread(fid,1,'uint16');
            EM_WaterColumn.SystemSerialNumber(ii)                = fread(fid,1,'uint16');
            EM_WaterColumn.NumberOfDatagrams(ii)                 = fread(fid,1,'uint16');
            EM_WaterColumn.DatagramNumbers(ii)                   = fread(fid,1,'uint16');
            EM_WaterColumn.NumberOfTransmitSectors(ii)           = fread(fid,1,'uint16'); %Ntx 
            EM_WaterColumn.TotalNumberOfReceiveBeams(ii)         = fread(fid,1,'uint16');
            EM_WaterColumn.NumberOfBeamsInThisDatagram(ii)       = fread(fid,1,'uint16'); %Nrx
            EM_WaterColumn.SoundSpeed(ii)                        = fread(fid,1,'uint16'); %SS 
            EM_WaterColumn.SamplingFrequency(ii)                 = fread(fid,1,'uint32'); %SF
            EM_WaterColumn.TXTimeHeave(ii)                       = fread(fid,1,'int16');
            EM_WaterColumn.TVGFunctionApplied(ii)                = fread(fid,1,'uint8'); %X
            EM_WaterColumn.TVGOffset(ii)                         = fread(fid,1,'int8'); %C
            EM_WaterColumn.ScanningInfo(ii)                      = fread(fid,1,'uint8');
            EM_WaterColumn.Spare1(ii)                            = fread(fid,1,'uint8');
            EM_WaterColumn.Spare2(ii)                            = fread(fid,1,'uint8');
            EM_WaterColumn.Spare3(ii)                            = fread(fid,1,'uint8');
            % repeat cycle #1: Ntx entries of 6 bits
                temp = ftell(fid);
                C = 6;
                Ntx = EM_WaterColumn.NumberOfTransmitSectors(ii);
                EM_WaterColumn.TiltAngle{ii}                     = fread(fid,Ntx,'int16',C-2);
                fseek(fid,temp+2,'bof'); % to next data type       
                EM_WaterColumn.CenterFrequency{ii}               = fread(fid,Ntx,'uint16',C-2);
                fseek(fid,temp+4,'bof'); % to next data type
                EM_WaterColumn.TransmitSectorNumber{ii}          = fread(fid,Ntx,'uint8',C-1);
                fseek(fid,temp+5,'bof'); % to next data type
                EM_WaterColumn.Spare{ii}                         = fread(fid,Ntx,'uint8',C-1);
                fseek(fid,1-C,'cof'); % we need to come back after last jump
            % repeat cycle #2: Nrx entries of a possibly variable number of bits. Using a for loop   
            Nrx = EM_WaterColumn.NumberOfBeamsInThisDatagram(ii);
            Ns = nan(1,Nrx);
            for jj=1:Nrx
                EM_WaterColumn.BeamPointingAngle{ii}(jj)             = fread(fid,1,'int16');
                EM_WaterColumn.StartRangeSampleNumber{ii}(jj)        = fread(fid,1,'uint16');
                EM_WaterColumn.NumberOfSamples{ii}(jj)               = fread(fid,1,'uint16'); %Ns
                EM_WaterColumn.DetectedRangeInSamples{ii}(jj)        = fread(fid,1,'uint16'); %DR
                EM_WaterColumn.TransmitSectorNumber2{ii}(jj)         = fread(fid,1,'uint8');
                EM_WaterColumn.BeamNumber{ii}(jj)                    = fread(fid,1,'uint8');
                Ns(jj) = EM_WaterColumn.NumberOfSamples{ii}(jj);
                EM_WaterColumn.SampleAmplitude{ii}{jj}               = fread(fid,Ns(jj),'int8');
            end
            % "spare byte if required to get even length (always 0 if used)"
            if floor((Nrx*10+sum(Ns))/2) == (Nrx*10+sum(Ns))/2
                % even so far, since ETX is 1 byte, add a spare here
                EM_WaterColumn.Spare4(ii)                            = fread(fid,1,'uint8');
            else
                % odd so far, since ETX is 1 bytes, no spare
                EM_WaterColumn.Spare4(ii) = NaN;
            end
            EM_WaterColumn.ETX(ii)                               = fread(fid,1,'uint8');
            EM_WaterColumn.CheckSum(ii)                          = fread(fid,1,'uint16');
            
            % ETX check
            if EM_WaterColumn.ETX(ii)~=3,
                error('wrong ETX value (EM_WaterColumn)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'WATER COLUMN DATAGRAM (107, 6BH)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end

        case 110 % NETWORK ATTITUDE VELOCITY DATAGRAM 110 (6EH)

            % datagrams counter
            if ~exist('EM_NetworkAttitude'), ii=1; else ii=size(EM_NetworkAttitude.TypeOfDatagram,2)+1; end

            % parsing
            EM_NetworkAttitude.NumberOfBytesInDatagram(ii)                    = nbDatag;
            EM_NetworkAttitude.STX(ii)                                        = stxDatag;
            EM_NetworkAttitude.TypeOfDatagram(ii)                             = typeDatag;
            EM_NetworkAttitude.EMModelNumber(ii)                              = emNumber;
            EM_NetworkAttitude.Date(ii)                                       = fread(fid,1,'uint32');
            EM_NetworkAttitude.TimeSinceMidnightInMilliseconds(ii)            = fread(fid,1,'uint32');
            EM_NetworkAttitude.NetworkAttitudeCounter(ii)                     = fread(fid,1,'uint16');
            EM_NetworkAttitude.SystemSerialNumber(ii)                         = fread(fid,1,'uint16');          
            EM_NetworkAttitude.NumberOfEntries(ii)                            = fread(fid,1,'uint16'); %N
            EM_NetworkAttitude.SensorSystemDescriptor(ii)                     = fread(fid,1,'int8'); 
            EM_NetworkAttitude.Spare(ii)                                      = fread(fid,1,'uint8'); 
            % repeat cycle: N entries of a variable number of bits. Using a for loop
            N = EM_NetworkAttitude.NumberOfEntries(ii);
            Nx = nan(1,N);
            for jj=1:N
                EM_NetworkAttitude.TimeInMillisecondsSinceRecordStart{ii}(jj)     = fread(fid,1,'uint16');
                EM_NetworkAttitude.Roll{ii}(jj)                                   = fread(fid,1,'int16');
                EM_NetworkAttitude.Pitch{ii}(jj)                                  = fread(fid,1,'int16');
                EM_NetworkAttitude.Heave{ii}(jj)                                  = fread(fid,1,'int16');
                EM_NetworkAttitude.Heading{ii}(jj)                                = fread(fid,1,'uint16');
                EM_NetworkAttitude.NumberOfBytesInInputDatagrams{ii}(jj)          = fread(fid,1,'uint8'); %Nx
                Nx(jj) = EM_NetworkAttitude.NumberOfBytesInInputDatagrams{ii}(jj);
                EM_NetworkAttitude.NetworkAttitudeInputDatagramAsReceived{ii}{jj} = fread(fid,Nx(jj),'uint8');
            end
            % "spare byte if required to get even length (always 0 if used)"
            if floor((N*11+sum(Nx))/2) == (N*11+sum(Nx))/2
                % even so far, since ETX is 1 byte, add a spare here
                EM_NetworkAttitude.Spare2(ii)                                    = fread(fid,1,'uint8');
            else
                % odd so far, since ETX is 1 bytes, no spare
                EM_NetworkAttitude.Spare2(ii) = NaN;
            end
            EM_NetworkAttitude.ETX(ii)                                           = fread(fid,1,'uint8');
            EM_NetworkAttitude.CheckSum(ii)                                      = fread(fid,1,'uint16');

            % ETX check
            if EM_NetworkAttitude.ETX(ii)~=3
                error('wrong ETX value (EM_NetworkAttitude)');
            end
            
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {'NETWORK ATTITUDE VELOCITY DATAGRAM 110 (110, 6EH)'};
                ALLfileinfo.parsedstatus(kk,1) = 1;
                ALLfileinfo.datagramcounter(kk,1) = ii;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end

        otherwise

            % typeDatag not recognized yet
           
            % file information
            if nargout
                kk = kk+1;
                ALLfileinfo.datagramtypenumber(kk,1) = typeDatag;
                ALLfileinfo.datagramtype(kk,1) = {sprintf('NOT SUPPORTED (%i, %sH)',typeDatag, dec2hex(typeDatag))};
                ALLfileinfo.parsedstatus(kk,1) = 0;
                ALLfileinfo.datagramcounter(kk,1) = NaN;
                ALLfileinfo.datagramsize(kk,1) = nbDatag;
                ALLfileinfo.synccounter(kk,1) = syncn;
            end
            
    end
            
    % reinitialize synccounter
    syncn = 0;
    
    % go to end of datagram
    fseek(fid,point+4+nbDatag,-1);     

end


%% saving data
fclose(fid);
save(MATfilename, '-regexp', 'EM\w*');
