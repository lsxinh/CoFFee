function CFF_multi_LOD_analysis(Z1_file,Z2_file,polygon,lod_method,U1_file,U2_file)
% CFF_multi_LOD_analysis(varargin)
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
%%%
% Alex Schimel, Deakin University
%%%


% read grid tif files:
[Z1,Z1_easting,Z1_northing] = CFF_readtif(Z1_file);
[Z2,Z2_easting,Z2_northing] = CFF_readtif(Z2_file);

% read uncertainty tif files:
[U1,U1_easting,U1_northing] = CFF_readtif(U1_file);
[U2,U2_easting,U2_northing] = CFF_readtif(U2_file);

% clip grids to polygons
[Z1,Z1_easting,Z1_northing] = CFF_clip_grid(Z1,Z1_easting,Z1_northing,polygon);
[Z2,Z2_easting,Z2_northing] = CFF_clip_grid(Z2,Z2_easting,Z2_northing,polygon);
[U1,U1_easting,U1_northing] = CFF_clip_grid(U1,U1_easting,U1_northing,polygon);
[U2,U2_easting,U2_northing] = CFF_clip_grid(U2,U2_easting,U2_northing,polygon);

% create dod from grids
[DOD,DOD_easting,DOD_northing,DOD_res] = CFF_create_DOD(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing,Z_res);

% now run multi-LOD analysis:
for t = 1:100
    
    lod_method = ['Ucrit',t];
    
    % create lod from std grids
    [LOD,LOD_easting,LOD_northing] = CFF_create_LOD(lod_method);
    
    % clip LOD
    [LOD,LOD_easting,LOD_northing] = CFF_clip_grid(LOD,LOD_easting,LOD_northing,polygon);
    
    [V_ero(t), V_dep(t)]  = CFF_LOD_analysis(DOD,LOD,interval);
    
end
