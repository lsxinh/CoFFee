function [Z,X,Y] = CFF_clip_grid(Z,X,Y,xv,yv)
% [Z,X,Y] = CFF_clip_grid(Z,X,Y,polygon)
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

% polygons in indices of grids:
xv = reshape(xv,1,numel(xv));
yv = reshape(yv,numel(yv),1);

% in easting northing coordinates:
poly = [easting(1,xv)', northing(yv,1)];

% build mask for each polygon
mask = nan(size(Z));
mask(inpolygon(easting,northing,poly(:,1),poly(:,2))) = 1;

% find the minimum and maximum of easting and northing

% compute and mask bathy diff
bathyDiff = bathyDiff .* bathyMask;




[Z1,Z1_easting,Z1_northing] = CFF_clip_grid(Z1,Z1_easting,Z1_northing,polygon);

