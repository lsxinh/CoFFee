function [M,easting,northing] = CFF_readtif(tiffile,tfwfile)
% CFF_new_function(varargin)
%
% DESCRIPTION
%
% read tif and tfw file
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
% ...
%
%%%
% Alex Schimel, Deakin University
%%%

% read tif file
M = imread(tiffile);
M = double(M);
M( M == min(M(:)) ) = NaN; % some tif have NaN values in the second layer. use imfinfo to figure out and

% read tfw file
tfw = csvread(tfwfile);

% turn tfw to easting/northing grid
row = [1:size(M,1)]';
col = 1:size(M,2);
[col,row] = meshgrid(col,row);
easting = tfw(1).*col + tfw(3).*row + tfw(5);
northing = tfw(2).*col + tfw(4).*row + tfw(6);

