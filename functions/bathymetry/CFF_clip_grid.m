function [Z2,X2,Y2] = CFF_clip_grid(Z,X,Y,xv,yv)
% [Z2,X2,Y2] = CFF_clip_grid(Z,X,Y,xv,yv)
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

% build mask
mask = nan(size(Z));
mask(inpolygon(X,Y,xv,yv)) = 1;

% apply mask
Z = Z.*mask;

% find limits of data:
dataZ = ~isnan(Z);

rows = double(any(dataZ,2));
irow_beg = find(rows,1,'first'); 
irow_end = find(rows,1,'last'); 

cols = double(any(dataZ,1));
icol_beg = find(cols,1,'first'); 
icol_end = find(cols,1,'last'); 

% output
Z2 = Z(irow_beg:irow_end,icol_beg:icol_end);
X2 = X(irow_beg:irow_end,icol_beg:icol_end);
Y2 = Y(irow_beg:irow_end,icol_beg:icol_end);


