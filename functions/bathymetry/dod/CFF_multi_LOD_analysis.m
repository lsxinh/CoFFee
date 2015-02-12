function CFF_multi_LOD_analysis(Z1_file,Z2_file,polygon,U1_file,U2_file)
% CFF_multi_LOD_analysis(Z1_file,Z2_file,polygon,U1_file,U2_file)
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


% read grid tif files:
[Z1,Z1_easting,Z1_northing] = CFF_readtif(Z1_file);
[Z2,Z2_easting,Z2_northing] = CFF_readtif(Z2_file);

% read uncertainty tif files:
[U1,U1_easting,U1_northing] = CFF_readtif(U1_file,Z1_file);
[U2,U2_easting,U2_northing] = CFF_readtif(U2_file,Z2_file);

% get polygons vertices:
xv = polygon(:,1);
yv = polygon(:,2);

% clip grids to polygons
[Z1,Z1_easting,Z1_northing] = CFF_clip_grid(Z1,Z1_easting,Z1_northing,xv,yv);
[Z2,Z2_easting,Z2_northing] = CFF_clip_grid(Z2,Z2_easting,Z2_northing,xv,yv);
[U1,U1_easting,U1_northing] = CFF_clip_grid(U1,U1_easting,U1_northing,xv,yv);
[U2,U2_easting,U2_northing] = CFF_clip_grid(U2,U2_easting,U2_northing,xv,yv);

% create dod from grids and propagated uncertainty grid
[DOD,DOD_easting,DOD_northing] = CFF_create_DOD(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing);
[DPU,DPU_easting,DPU_northing] = CFF_create_DPU(U1,U1_easting,U1_northing,U2,U2_easting,U2_northing);

% grids should be co-registered, but just in case of:
[DOD,DPU,X,Y] = CFF_coregister_grids(DOD,DOD_easting,DOD_northing,DPU,DPU_easting,DPU_northing);

% now run multi-LOD analysis on several levels of confidence intervals
for i = 1:100
    
    % LOD analysis
    [v_bud(i),v_ero(i),v_dep(i),a_ero(i),a_dep(i),us_v_ero(i),us_v_dep(i),up_v_ero(i),up_v_dep(i)] = CFF_LOD_analysis(DOD,DPU,X,Y,i);
    
end

% display
figure;

plot([1:100],v_dep,'Color',[0.4 0.4 0.4],'LineWidth',2)
hold on
plot([1:100],v_bud,'Color',[0 0 0],'LineWidth',2)
plot([1:100],v_ero,'Color',[0.7 0.7 0.7],'LineWidth',2)
legend('deposition','net','erosion')

plot([1:100],v_ero-us_v_ero,'--','Color',[0.7 0.7 0.7],'LineWidth',2)
plot([1:100],v_dep+us_v_dep,'--','Color',[0.4 0.4 0.4],'LineWidth',2)
us_v_ero(v_ero+us_v_ero > 0) = NaN;
plot([1:100],v_ero+us_v_ero,'--','Color',[0.7 0.7 0.7],'LineWidth',2)
us_v_dep(v_dep-us_v_dep < 0) = NaN;
plot([1:100],v_dep-us_v_dep,'--','Color',[0.4 0.4 0.4],'LineWidth',2)

vmax = max([max(abs(v_ero-us_v_ero)), max(abs(v_dep+us_v_dep))]);
stem([0.14,0.14],[-vmax,vmax],'k.','LineWidth',2)

grid on

xlabel('[1:100] (m)')
ylabel('erosion (m^3)                     deposition (m^3)')
title(['Volumes above threshold for ' bathyDiffLabel ' on ' bathyMaskLabel ' area'])
ylim([-vmax vmax])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0.25 0.25 30 20]);
CFF_nice_easting_northing(5)
%print('-dpng','-r600',[bathyDiffLabel '_' bathyMaskLabel '_volumeAboveThreshold.png'])



        
% bathy intervals to compute interval volumes
interval = [-1:0.01:1]; %m

