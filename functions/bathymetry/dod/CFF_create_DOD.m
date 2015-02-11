function [DOD,Xout,Yout] = CFF_create_DOD(Z1,X1,Y1,Z2,X2,Y2,varargin)
% [DOD,Xout,Yout] = CFF_create_DOD(Z1,X1,Y1,Z2,X2,Y2,varargin)
%
% DESCRIPTION
%
% Compute the difference between the two grids. Includes a call to
% co-register the grids on the same X,Y vectors in case it's not the case.
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

% coregister the grids depending on varargin
if nargin==6
    % unspecified
    [Z1out,Z2out,Xout,Yout] = CFF_coregister_grids(Z1,X1,Y1,Z2,X2,Y2);
elseif nargin==7
    % one value specified, pass it on
    [Z1out,Z2out,Xout,Yout] = CFF_coregister_grids(Z1,X1,Y1,Z2,X2,Y2,varargin{1});
elseif nargin==8
    % two values specified, pass them on
    [Z1out,Z2out,Xout,Yout] = CFF_coregister_grids(Z1,X1,Y1,Z2,X2,Y2,varargin{1},varargin{2});
end

% compute difference of DEM
DOD = Z2out-Z1out;



