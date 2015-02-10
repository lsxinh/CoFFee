function [LOD,LOD_easting,LOD_northing] = CFF_create_LOD(lod_method)
% [LOD,LOD_easting,LOD_northing] = CFF_create_LOD(lod_method)
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
% std1 = '.\DATA\WH1_uncertaintyZ_50cm_UTM54S_p.tif';
% std2 = '.\DATA\WH2_uncertaintyZ_50cm_UTM54S_p.tif';
%
%%%
% Alex Schimel, Deakin University
%%%


% read std tif if they exist
if nargin>3
    
    std1 = varargin{1};
    std2 = varargin{2};
    
    % getting tfw names
    [pathstr1,name1,ext1] = fileparts(std1);
    stw1 = [pathstr1 filesep name1 '.tfw'];
    [pathstr2,name2,ext2] = fileparts(std2);
    stw2 = [pathstr2 filesep name2 '.tfw'];
    
    % read tif and tfw
    if exist(stw1,'file')
        [ST1,ST1_easting,ST1_northing] = CFF_readtif(std1,stw1);
    else
        [ST1,ST1_easting,ST1_northing] = CFF_readtif(std1);
    end
    if exist(stw2,'file')
        [ST2,ST2_easting,ST2_northing] = CFF_readtif(std2,stw2);
    else
        [ST2,ST2_easting,ST2_northing] = CFF_readtif(std2);
    end
    
else
    
    % if we don't have some uncertainty value, use a fixed value, taken
    % from ref area: eg 0.14. So test
    
end

% is ST1 has no georeference but is the same size as Z1, use Z1 georeference
if ST1_easting(1)==1 & Z1_easting(1)~=1 & all(size(ST1)==size(Z1))
    ST1_easting = Z1_easting;
    ST1_northing = Z1_northing;
end
% is ST2 has no georeference but is the same size as Z2, use Z2 georeference
if ST2_easting(1)==1 & Z2_easting(1)~=1 & all(size(ST2)==size(Z2))
    ST2_easting = Z2_easting;
    ST2_northing = Z2_northing;
end
    

ST1res = ST1_easting(1,2)-ST1_easting(1,1);
ST2res = ST2_easting(1,2)-ST2_easting(1,1);

% get coordinates for the max extent
min_easting = min([Z1_easting(:);Z2_easting(:);ST1_easting(:);ST2_easting(:)]);
max_easting = max([Z1_easting(:);Z2_easting(:);ST1_easting(:);ST2_easting(:)]);
min_northing = min([Z1_northing(:);Z2_northing(:);ST1_northing(:);ST2_northing(:)]);
max_northing = max([Z1_northing(:);Z2_northing(:);ST1_northing(:);ST2_northing(:)]);
[easting,northing] = meshgrid([min_easting:Zres:max_easting],[max_northing:-Zres:min_northing]);



firstcol = find(easting(1,:)==ST1_easting(1));
firstrow = find(northing(:,1)==ST1_northing(1));
Zblank = nan(size(easting));
Zblank( firstrow:firstrow+size(ST1,1)-1 , firstcol:firstcol+size(ST1,2)-1 ) = ST1;
bigST1 = Zblank;

firstcol = find(easting(1,:)==ST2_easting(1));
firstrow = find(northing(:,1)==ST2_northing(1));
Zblank = nan(size(easting));
Zblank( firstrow:firstrow+size(ST2,1)-1 , firstcol:firstcol+size(ST2,2)-1 ) = ST2;
bigST2 = Zblank;


% compute LOD
bigSig = sqrt( bigST1.^2 + bigST2.^2 );

t = 1; % 1 sigma, 68% confidence limit
t = 1.96; % 2 sigma, 95% confidence limit
Ucrit = t.*bigSig; %threshold to use, as per Milan et al., 2011; Brasington et al., 2003; Lane et al., 2003; lallias tacon 2013; Wheaton et al., 2010


%"minimum level of detection threshold (LOD)" lOD = Ucrit
%"LoD grid" (Carley et al., 2012)




% using 1.96 as Ucrit threshold?



% output co-registered grids, DoD, LoD


