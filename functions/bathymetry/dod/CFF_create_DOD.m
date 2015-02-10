function [DOD,DOD_easting,DOD_northing,DOD_res] = CFF_create_DOD(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing,Z_res)
% [DOD,DOD_easting,DOD_northing,DOD_res] = CFF_create_DOD(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing,Z_res)
%
% DESCRIPTION
%
% use as template for a new function
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
% - varagin
%
% OUTPUT VARIABLES
%
% - NA
%
% RESEARCH NOTES
%
% ...
%
% NEW FEATURES
%
% YYYY-MM-DD: second version.
% YYYY-MM-DD: first version.
%
% EXAMPLE
%
% tif1 = '.\DATA\WH1_Z_50cm_UTM54S_LAT_p1.tif';
% tif2 = '.\DATA\WH2_Z_50cm_UTM54S_LAT_p1.tif';

%
%%%
% Alex Schimel, Deakin University
%%%

% getting tfw names
[pathstr1,name1,ext1] = fileparts(tif1);
tfw1 = [pathstr1 filesep name1 '.tfw'];

[pathstr2,name2,ext2] = fileparts(tif2);
tfw2 = [pathstr2 filesep name2 '.tfw'];

% read tif and tfw
if exist(tfw1,'file')
    [Z1,Z1_easting,Z1_northing] = CFF_readtif(tif1,tfw1);
else
    [Z1,Z1_easting,Z1_northing] = CFF_readtif(tif1);
end
if exist(tfw2,'file')
    [Z2,Z2_easting,Z2_northing] = CFF_readtif(tif2,tfw2);
else
    [Z2,Z2_easting,Z2_northing] = CFF_readtif(tif2);
end



% of now some files may have proper easting/northing, some may not. deal
% with it by testing for size and assuming they are the same.

% if Z1 has no georeference but is the same size as Z2, use Z2 georeference
if Z1_easting(1)==1 & Z2_easting(1)~=1 & all(size(Z1)==size(Z2))
    Z1_easting = Z2_easting;
    Z1_northing = Z2_northing;
end
% if Z2 has no georeference but is the same size as Z1, use Z1 georeference
if Z1_easting(1)~=1 & Z2_easting(1)==1 & all(size(Z1)==size(Z2))
    Z2_easting = Z1_easting;
    Z2_northing = Z1_northing;
end

% get resolutions:
Z1res = Z1_easting(1,2)-Z1_easting(1,1);
Z2res = Z2_easting(1,2)-Z2_easting(1,1);


% at this stage, I can only deal with all resolutions being the same
if all(Z1res == [Z2res,ST1res,ST2res])
    Zres = Z1res;
else
    error;
end

% now work on putting everything on the same exact grid:

% get coordinates for the max extent
min_easting = min([Z1_easting(:);Z2_easting(:)]);
max_easting = max([Z1_easting(:);Z2_easting(:)]);
min_northing = min([Z1_northing(:);Z2_northing(:)]);
max_northing = max([Z1_northing(:);Z2_northing(:)]);
[easting,northing] = meshgrid([min_easting:Zres:max_easting],[max_northing:-Zres:min_northing]);

% change grids to the max extent
firstcol = find(easting(1,:)==Z1_easting(1));
firstrow = find(northing(:,1)==Z1_northing(1));
Zblank = nan(size(easting));
Zblank( firstrow:firstrow+size(Z1,1)-1 , firstcol:firstcol+size(Z1,2)-1 ) = Z1;
bigZ1 = Zblank;

firstcol = find(easting(1,:)==Z2_easting(1));
firstrow = find(northing(:,1)==Z2_northing(1));
Zblank = nan(size(easting));
Zblank( firstrow:firstrow+size(Z2,1)-1 , firstcol:firstcol+size(Z2,2)-1 ) = Z2;
bigZ2 = Zblank;


% stacking all bathy to measure for each cell: count of bathy dataset with
% value, depth range, and whether all data are available
temp = cat(3,bigZ1,bigZ2);
bathyRange = max(temp,[],3) - min(temp,[],3);
bathyCount = sum(double(~isnan(temp)),3);
bathyAll = double(bathyCount==2);

% compute difference
DoD = bigZ2-bigZ1;



