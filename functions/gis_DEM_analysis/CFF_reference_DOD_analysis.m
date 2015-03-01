function [DODmean,DODstd,DOD,X,Y] = CFF_reference_DOD_analysis(Z1_file,Z2_file,polygon)
% [DODmean,DODstd,DOD,X,Y] = CFF_reference_DOD_analysis(Z1_file,Z2_file,polygon)
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

% get polygons vertices:
xv = polygon(:,1);
yv = polygon(:,2);

% read grid tif files:
[Z1,Z1_easting,Z1_northing] = CFF_read_tif(Z1_file);
[Z2,Z2_easting,Z2_northing] = CFF_read_tif(Z2_file);

% clip grids to polygon
[Z1,Z1_easting,Z1_northing] = CFF_clip_raster(Z1,Z1_easting,Z1_northing,xv,yv);
[Z2,Z2_easting,Z2_northing] = CFF_clip_raster(Z2,Z2_easting,Z2_northing,xv,yv);

% coregister grids
[Z1,Z2,X,Y] = CFF_coregister_rasters(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing);

% create dod from grids 
DOD = CFF_calculate_DOD(Z1,Z2);

% get mean and standard deviation of DOD over reference area
[DODmean,DODstd] = CFF_nanstat3(DOD(:),1);

