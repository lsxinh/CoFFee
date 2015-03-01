function CFF_single_LOD_analysis(Z1_file,Z2_file,polygon,varargin)
% CFF_single_LOD_analysis(Z1_file,Z2_file,polygon,varargin)
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

% get polygons vertices:
xv = polygon(:,1);
yv = polygon(:,2);

tic
% read grid tif files:
[Z1,Z1_easting,Z1_northing] = CFF_read_tif(Z1_file);
[Z2,Z2_easting,Z2_northing] = CFF_read_tif(Z2_file);

toc
tic

% clip grids to polygon
[Z1,Z1_easting,Z1_northing] = CFF_clip_raster(Z1,Z1_easting,Z1_northing,xv,yv);
[Z2,Z2_easting,Z2_northing] = CFF_clip_raster(Z2,Z2_easting,Z2_northing,xv,yv);

toc

% create dod from grids
[DOD,DOD_easting,DOD_northing] = CFF_calculate_DOD(Z1,Z1_easting,Z1_northing,Z2,Z2_easting,Z2_northing);

if nargin>4
    
    % we have uncertainty files
    U1_file = varargin{1};
    U2_file = varargin{2};
    
    % read uncertainty tif files:
    [U1,U1_easting,U1_northing] = CFF_read_tif(U1_file,Z1_file);
    [U2,U2_easting,U2_northing] = CFF_read_tif(U2_file,Z2_file);
    
    % clip uncertainty grids to polygon
    [U1,U1_easting,U1_northing] = CFF_clip_raster(U1,U1_easting,U1_northing,xv,yv);
    [U2,U2_easting,U2_northing] = CFF_clip_raster(U2,U2_easting,U2_northing,xv,yv);
    
    % create propagated uncertainty grid
    [DPU,DPU_easting,DPU_northing] = CFF_calculate_DPU(U1,U1_easting,U1_northing,U2,U2_easting,U2_northing);
    
    % DOD and DPU should be co-registered, but just in case of:
    [DOD,DPU,DOD_easting,DOD_northing] = CFF_coregister_rasters(DOD,DOD_easting,DOD_northing,DPU,DPU_easting,DPU_northing);
    
    % now run multi-LOD analysis
    sigma = [0:1:30];
    for i = 1:length(sigma)
        
        % we threshold at a spatially variable value, a factor of the
        % grid std, propagated in quadrature from the individual grids.
        uncertainty = DPU;
        threshold = sigma(i).*DPU;
        [v_bud(i),v_ero(i),v_dep(i),a_ero(i),a_dep(i),us_v_ero(i),us_v_dep(i),up_v_ero(i),up_v_dep(i)] = CFF_LOD_volumes(DOD,DOD_easting,DOD_northing,threshold,uncertainty);
        
    end
    
elseif nargin >3
    
    % we have a single uncertainty value
    uncertainty = varargin{1};
    
    % now run multi-LOD analysis
    % we threshold at a constant value, a factor of some std (reference area)
    sigma = [0:0.1:5];
    threshold = sigma.*uncertainty;
    for i = 1:length(threshold)
        [v_bud(i),v_ero(i),v_dep(i),a_ero(i),a_dep(i),us_v_ero(i),us_v_dep(i),up_v_ero(i),up_v_dep(i)] = CFF_LOD_volumes(DOD,DOD_easting,DOD_northing,threshold(i),uncertainty);
    end
    
end



% display
figure;

plot(threshold,v_dep,'Color',[0.4 0.4 0.4],'LineWidth',2)
hold on
plot(threshold,v_bud,'Color',[0 0 0],'LineWidth',2)
plot(threshold,v_ero,'Color',[0.7 0.7 0.7],'LineWidth',2)
legend('deposition','net','erosion')

% erosion uncertainty
plot(threshold,v_ero-us_v_ero,'--','Color',[0.7 0.7 0.7],'LineWidth',2)
us_v_ero(v_ero+us_v_ero > 0) = NaN;
plot(threshold,v_ero+us_v_ero,'--','Color',[0.7 0.7 0.7],'LineWidth',2)

% deposition uncertainty
plot(threshold,v_dep+us_v_dep,'--','Color',[0.4 0.4 0.4],'LineWidth',2)
us_v_dep(v_dep-us_v_dep < 0) = NaN;
plot(threshold,v_dep-us_v_dep,'--','Color',[0.4 0.4 0.4],'LineWidth',2)

% value at uncertainty level
vmax = max([max(abs(v_ero-us_v_ero)), max(abs(v_dep+us_v_dep))]);
stem([uncertainty,uncertainty],[-vmax,vmax],'k.','LineWidth',2)

grid on
xlabel('threshold (m)')
ylabel('erosion (m^3)                     deposition (m^3)')
title(['Volumes above threshold'])
ylim([-vmax vmax])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0.25 0.25 30 20]);
CFF_nice_easting_northing(5)
print('-dpng','-r600','volumeAboveThreshold.png')


