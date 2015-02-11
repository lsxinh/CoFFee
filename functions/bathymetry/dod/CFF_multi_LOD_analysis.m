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
[U1,U1_easting,U1_northing] = CFF_readtif(U1_file,Z1_file);
[U2,U2_easting,U2_northing] = CFF_readtif(U2_file,Z2_file);

% get polygons vertices:
xv = polygon(:,1);
yv = polygon(:,2);

% clip grids to polygons
[Z1,Z1_easting,Z1_northing] = CFF_clip_grid(Z1,Z1_easting,Z1_northing,xv,yv);
[Z2,Z2_easting,Z2_northing] = CFF_clip_grid(Z2,Z2_easting,Z2_northing,xv,yv);
[U1,U1_easting,U1_northing] = CFF_clip_grid(U1,U1_easting,U1_northing,xv,yv);
[U2,U2_easting,U2_northing] = CFF_clip_grid(U2,U2_easting,U2_northing,xv,yv);

% create dod from grids
[DOD,DOD_easting,DOD_northing] = CFF_create_DOD(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing);

% now run multi-LOD analysis:
for t = 1:100

    % create lod from std grids
    [LOD,LOD_easting,LOD_northing] = CFF_create_LOD('Ucrit',t,U1,U1_easting,U1_northing,U2,U2_easting,U2_northing);
    
    % grids should be co-registered, but just in case of:
    [DOD,LOD,X,Y] = CFF_coregister_grids(DOD,DOD_easting,DOD_northing,LOD,LOD_easting,LOD_northing);
    
    % LOD analysis
    
    % bathy intervals to compute interval volumes
    interval = [-1:0.01:1]; %m
    [V_ero(t), V_dep(t)]  = CFF_LOD_analysis(DOD,LOD,X,Y,interval);
    
end
