function FPBS = CFF_concatenate_FPBS(FPBS1,FPBS2)
% FPBS = CFF_concatenate_FPBS(FPBS1,FPBS2)
%
% DESCRIPTION
%
% This function merges two sets of multibeam data in the FPBS format into
% one.
%
% USE
%
% ...
%
% PROCESSING SUMMARY
% 
% - ...
% - ...
% - ...
%
% INPUT VARIABLES
%
% - FPBS1
% - FPBS2
%
% OUTPUT VARIABLES
%
% - FPBS
%
% RESEARCH NOTES
%
% ...
%
% NEW FEATURES
%
% 2013-09-30: first version.
%
% EXAMPLE
%
% ...
%
%%%
% Alex Schimel, Deakin University
%%%


% save the first FPBS as is
FPBS = FPBS1;

% get number of entries in each table of that first FPBS (should be equal
% to last ID value)
nf1 = length(FPBS.File.ID);
np1 = length(FPBS.Ping.ID);
nb1 = length(FPBS.Beam.ID);
ns1 = length(FPBS.Samp.ID);

% udpate ID in FPBS2 (just stack numbers)
FPBS2.File.ID = FPBS2.File.ID + nf1;
FPBS2.Ping.ID = FPBS2.Ping.ID + np1;
FPBS2.Beam.ID = FPBS2.Beam.ID + nb1;
FPBS2.Samp.ID = FPBS2.Samp.ID + ns1;

% udpate Index in FPBS2 (add number of rows from previous table)
FPBS2.Ping.Index = FPBS2.Ping.Index + nf1;
FPBS2.Beam.Index = FPBS2.Beam.Index + np1;
FPBS2.Samp.Index = FPBS2.Samp.Index + nb1;

% retrieve the structure field names
fileFields = fieldnames(FPBS.File);
pingFields = fieldnames(FPBS.Ping);
beamFields = fieldnames(FPBS.Beam);
sampFields = fieldnames(FPBS.Samp);

% Loop over the field-names and apply the function to each field
for n = 1:length(fileFields)
    FPBS.File.(fileFields{n}) = [ FPBS.File.(fileFields{n}) , FPBS2.File.(fileFields{n}) ];
end
for n = 1:length(pingFields)
    FPBS.Ping.(pingFields{n}) = [ FPBS.Ping.(pingFields{n}) , FPBS2.Ping.(pingFields{n}) ];
end
for n = 1:length(beamFields)
    FPBS.Beam.(beamFields{n}) = [ FPBS.Beam.(beamFields{n}) , FPBS2.Beam.(beamFields{n}) ];
end
for n = 1:length(sampFields)
    FPBS.Samp.(sampFields{n}) = [ FPBS.Samp.(sampFields{n}) , FPBS2.Samp.(sampFields{n}) ];
end

