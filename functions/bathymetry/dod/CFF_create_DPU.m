function [DPU,X,Y] = CFF_create_DPU(U1,X1,Y1,U2,X2,Y2,varargin)
% [DPU,X,Y] = CFF_create_DPU(U1,X1,Y1,U2,X2,Y2,varargin)
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
% std1 = '.\DATA\WH1_uncertaintyZ_50cm_UTM54S_p.tif';
% std2 = '.\DATA\WH2_uncertaintyZ_50cm_UTM54S_p.tif';
%
%%%
% Alex Schimel, Deakin University
%%%

% coregister the grids depending on varargin
if nargin==6
    % unspecified
    [U1,U2,X,Y] = CFF_coregister_grids(U1,X1,Y1,U2,X2,Y2);
elseif nargin==7
    % one value specified, pass it on
    [U1,U2,X,Y] = CFF_coregister_grids(U1,X1,Y1,U2,X2,Y2,varargin{1});
elseif nargin==8
    % two values specified, pass them on
    [U1,U2,X,Y] = CFF_coregister_grids(U1,X1,Y1,U2,X2,Y2,varargin{1},varargin{2});
end

% Uncertainty propagated in quadrature:
DPU = sqrt( U1.^2 + U2.^2 );



