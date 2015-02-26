function [DODmean,DODstd,DOD,DOD_easting,DOD_northing] = CFF_reference_DOD_analysis(Z1_file,Z2_file,polygon)
% [DODmean,DODstd,DOD,DOD_easting,DOD_northing] = CFF_reference_DOD_analysis(Z1_file,Z2_file,polygon)
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
[Z1,Z1_easting,Z1_northing] = CFF_readtif(Z1_file);
[Z2,Z2_easting,Z2_northing] = CFF_readtif(Z2_file);

% clip grids to polygon
[Z1,Z1_easting,Z1_northing] = CFF_clip_grid(Z1,Z1_easting,Z1_northing,xv,yv);
[Z2,Z2_easting,Z2_northing] = CFF_clip_grid(Z2,Z2_easting,Z2_northing,xv,yv);

% create dod from grids 
[DOD,DOD_easting,DOD_northing] = CFF_create_DOD(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing);

% get mean and standard deviation of DOD over reference area
[DODmean,DODstd] = CFF_nanstat3(DOD(:),1);