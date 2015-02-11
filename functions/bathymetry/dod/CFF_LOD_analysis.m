function [V_ero, V_dep] = CFF_LOD_analysis(DOD,LOD,X,Y,interval)
% [V_ero, V_dep] = CFF_LOD_analysis(DOD,LOD,interval)
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

% cell resolutions in X and Y and area
Xres = X(1,2)-X(1,1);
Yres = Y(1,1)-Y(2,1);
cellArea = Xres.*Yres;

% erosion
DOD_ero_mask = DOD<-LOD;
volumeErodedAboveLOD  = CFF_nansum3(CFF_nansum3(DOD .* DOD_ero_mask)).*cellArea;
areaErodedAboveLOD    = sum(sum(double(DOD_ero_mask))).*cellArea;

uncertaintyOfVolumeErodedAboveLOD = CFF_nansum3(CFF_nansum3(LOD .* DOD_ero_mask)).*cellArea; % sum of LODs for cells above LOD
uncertaintyOfVolumeErodedAboveLOD = sqrt( sum(sum(double(DOD_ero_mask))) .* betweenGridsPrecision.^2 ).*cellArea;

% deposition
DOD_dep_mask = DOD>LOD;
volumeDepositedAboveLOD  = CFF_nansum3(CFF_nansum3(DOD .* DOD_dep_mask)).*cellArea;
areaDepositedAboveLOD   = sum(sum(double(DOD_dep_mask))).*cellArea;

% budget
volumeBudgetAboveLOD = volumeDepositedAboveLOD + volumeErodedAboveLOD;

% uncertainty


uncertaintyOfVolumeDepositedAboveLOD = betweenGridsPrecision.*areaDepositedAboveLOD;
uncertaintyOfVolumeDepositedAboveLOD = sqrt( sum(sum(double(DOD_dep_mask))) .* betweenGridsPrecision.^2 ).*cellArea;



