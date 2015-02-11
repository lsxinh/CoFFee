function [M,easting,northing] = CFF_readtif(tif_file,varargin)
% CFF_new_function(varargin)
%
% DESCRIPTION
%
% read tif and tfw file. If tfw file unavailable, returns row and col
% number.
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
% - tif_file
% - varagin{1}: tfw file name
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
% 2015-02-10: tfw now optional.
% YYYY-MM-DD: first version.
%
% EXAMPLE
%
% ...
%
%%%
% Alex Schimel, Deakin University
%%%

% read tif file
M = imread(tif_file);
M = double(M);

% use minimum value in tiff as NaN value for now
% some tif have NaN values in the second layer. use imfinfo to figure out
% and change this code
M( M == min(M(:)) ) = NaN;

% indices grid
row = [1:size(M,1)]';
col = 1:size(M,2);
[col,row] = meshgrid(col,row);

% now find the tfw file
if nargin>1
    % check extension
    [pathstr,name,ext] = fileparts(varargin{1});
    if strcmp(ext,'.tfw')
        % input is a tfw file
        tfw_file = varargin{1};
    else
        % input is not a tfw file. Could be a tif file, which associated
        % tfw we want. Change extension
        tfw_file = [pathstr filesep name '.tfw'];
    end
else
    % try find a tfw associated with input tif
    [pathstr,name,ext] = fileparts(tif_file);
    tfw_file = [pathstr filesep name '.tfw'];
end

% now check existence of tfw file and read it
if exist(tfw_file,'file')
    % read
    tfw = csvread(tfw_file);
    % turn tfw to easting/northing grid
    easting = tfw(1).*col + tfw(3).*row + tfw(5);
    northing = tfw(2).*col + tfw(4).*row + tfw(6);
else
    % if tfw is not available, we'll use row and col for easting and
    % northing by default
    warning('Could not find a .tfw file. Exporting grid indices as easting/northing')
    easting = col;
    northing = flipud(row);
end

