function [v_bud,v_ero,v_dep,a_ero,a_dep,us_v_ero,us_v_dep,up_v_ero,up_v_dep] = CFF_LOD_volumes(DOD,X,Y,LOD,UNC)
% [v_bud,v_ero,v_dep,a_ero,a_dep,us_v_ero,us_v_dep,up_v_ero,up_v_dep] = CFF_LOD_volumes(DOD,X,Y,LOD,UNC)
%
% DESCRIPTION
%
% Calculate volumes eroded and deposited from two DEMS.
% using two co-registered grids: a difference of DEMs (DOD) and a threshold
% grid (LOD). ALso use uncertainty to output intervals of confidence.
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

% cell resolutions in X and Y and area
Xres = X(1,2)-X(1,1);
Yres = Y(1,1)-Y(2,1);
cellArea = Xres.*Yres;

%% erosion
DOD_ero_mask = DOD<-LOD;

% volume
v_ero  = CFF_nansum3(CFF_nansum3(DOD .* DOD_ero_mask .* cellArea));

% area
a_ero  = sum(sum(double(DOD_ero_mask))).*cellArea;

% The volume uncertainty of each cell is given by the product of the cell
% area by the cell uncertainty.
u_v_ero = UNC .* DOD_ero_mask .* cellArea;

% uncertainty in natural sum:
us_v_ero = CFF_nansum3(CFF_nansum3(u_v_ero)); 

% propagated uncertainty:
up_v_ero = sqrt(CFF_nansum3(CFF_nansum3(u_v_ero.^2)));

%% deposition
DOD_dep_mask = DOD>LOD;

% volume
v_dep  = CFF_nansum3(CFF_nansum3(DOD .* DOD_dep_mask .* cellArea));

% area
a_dep   = sum(sum(double(DOD_dep_mask))).*cellArea;

% The volume uncertainty of a cell is given by the product of the cell area
% by the cell DPU.
u_v_dep = UNC .* DOD_dep_mask .* cellArea;

% uncertainty in natural sum:
us_v_dep = CFF_nansum3(CFF_nansum3(u_v_dep)); 

% propagated uncertainty:
up_v_dep = sqrt(CFF_nansum3(CFF_nansum3(u_v_dep.^2)));

%% budget
v_bud = v_dep + v_ero;
