function [M,easting,northing] = CFF_readtif(tiffile,varargin)
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
% - tiffile
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
M = imread(tiffile);
M = double(M);
M( M == min(M(:)) ) = NaN; % some tif have NaN values in the second layer. use imfinfo to figure out and

% getting dimensions
row = [1:size(M,1)]';
col = 1:size(M,2);
[col,row] = meshgrid(col,row);

% read tfw file if it exists and get easting/northing
if nargin>1
    
    % read by replacing ext with tfw (to allow inputing other tif files to grab geo info from)
    [pathstr,name,ext] = fileparts(varargin{1});
    tfw = csvread([pathstr filesep name '.tfw']);
    
    % turn tfw to easting/northing grid
    easting = tfw(1).*col + tfw(3).*row + tfw(5);
    northing = tfw(2).*col + tfw(4).*row + tfw(6);
    
else
    
    % try find a tfw
    [pathstr,name,ext] = fileparts(tiffile);
    
    try
        
        tfw = csvread([pathstr filesep name '.tfw']);
        
        % turn tfw to easting/northing grid
        easting = tfw(1).*col + tfw(3).*row + tfw(5);
        northing = tfw(2).*col + tfw(4).*row + tfw(6);
        
    catch
        
        % if tfw is not available, start from 1
        easting = col;
        northing = flipud(row);
        
    end
    
end



