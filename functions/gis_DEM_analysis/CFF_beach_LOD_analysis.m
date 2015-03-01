function CFF_beach_LOD_analysis(Z1_file,Z2_file,reference_polygon,main_polygon,analysis_polygons)


%% read grid tif files:
[Z1,Z1_easting,Z1_northing] = CFF_read_tif(Z1_file);
[Z2,Z2_easting,Z2_northing] = CFF_read_tif(Z2_file);


%% get uncertainty from reference area

% % get polygons vertices:
% xv = reference_polygon(:,1);
% yv = reference_polygon(:,2);
% 
% % clip grids to polygon
% [ref_Z1,ref_X1,ref_Y1] = CFF_clip_raster(Z1,Z1_easting,Z1_northing,xv,yv);
% [ref_Z2,ref_X2,ref_Y2] = CFF_clip_raster(Z2,Z2_easting,Z2_northing,xv,yv);
% 
% % coregister grids
% [ref_Z1,ref_Z2,ref_X,ref_Y] = CFF_coregister_rasters(ref_Z1,ref_X1,ref_Y1,ref_Z2,ref_X2,ref_Y2);
% 
% % create dod from grids 
% ref_DOD = CFF_calculate_DOD(ref_Z1,ref_Z2);
% 
% % get mean and standard deviation of DOD over reference area
% [temp,uncertainty] = CFF_nanstat3(ref_DOD(:),1);

load uncertainty


%% clip to main polygon

% get polygons vertices:
xv = main_polygon(:,1);
yv = main_polygon(:,2);

% clip grids to polygon
[Z1,X1,Y1] = CFF_clip_raster(Z1,Z1_easting,Z1_northing,xv,yv);
[Z2,X2,Y2] = CFF_clip_raster(Z2,Z2_easting,Z2_northing,xv,yv);

% coregister grids
[Z1,Z2,X,Y] = CFF_coregister_rasters(Z1,X1,Y1,Z2,X2,Y2);


%% multi LOD analysis on main polygon

% create dod from grids 
DOD = CFF_calculate_DOD(Z1,Z2);

% now run multi-LOD analysis
% we threshold at a constant value, a factor of some std (reference area)
main_sigmafactor = [0:0.1:5];
main_threshold = main_sigmafactor.*uncertainty;
for i = 1:length(main_threshold)
    main_volumes(i) = CFF_LOD_volumes(DOD,X,Y,main_threshold(i),uncertainty);
end
    
save main_volumes main_volumes
    
%% single LOD analysis on multiple polygons

for ii = 1:length(analysis_polygons)
    
    ii
    
    % get polygon vertices:
    xv = analysis_polygons{ii}(:,1);
    yv = analysis_polygons{ii}(:,2);

    % clip grids to polygon
    [poly_Z1,poly_X1,poly_Y1] = CFF_clip_raster(Z1,X,Y,xv,yv);
    [poly_Z2,poly_X2,poly_Y2] = CFF_clip_raster(Z2,X,Y,xv,yv);
    
     % coregister grids
    [poly_Z1,poly_Z2,poly_X,poly_Y] = CFF_coregister_rasters(poly_Z1,poly_X1,poly_Y1,poly_Z2,poly_X2,poly_Y2);
    
    % create dod from grids
    poly_DOD = CFF_calculate_DOD(poly_Z1,poly_Z2);

    % now run multi-LOD analysis
    volumes_0(ii) = CFF_LOD_volumes(poly_DOD,poly_X,poly_Y,0,uncertainty);
    volumes_1(ii) = CFF_LOD_volumes(poly_DOD,poly_X,poly_Y,1,uncertainty);
    volumes_196(ii) = CFF_LOD_volumes(poly_DOD,poly_X,poly_Y,1.96,uncertainty);

end

save volumes volumes_0 volumes_1 volumes_196







% display
volumeNetChange = [main_volumes(:).volumeNetChange];
volumeEroded = [main_volumes(:).volumeEroded];
volumeDeposited = [main_volumes(:).volumeDeposited];
uncertaintyVolumeEroded_sum = [main_volumes(:).uncertaintyVolumeEroded_sum];
uncertaintyVolumeDeposited_sum = [main_volumes(:).uncertaintyVolumeDeposited_sum];
uncertaintyVolumeEroded_propagated = [main_volumes(:).uncertaintyVolumeEroded_propagated];
uncertaintyVolumeDeposited_propagated = [main_volumes(:).uncertaintyVolumeDeposited_propagated];
areaEroded = [main_volumes(:).areaEroded];
areaDeposited = [main_volumes(:).areaDeposited];
areaTotalChange = [main_volumes(:).areaTotalChange];
areaTotal = [main_volumes(:).areaTotal];

figure;

plot(main_threshold, volumeDeposited, 'Color',[0.4 0.4 0.4],'LineWidth',2)
hold on
plot(main_threshold, volumeNetChange, 'Color',[0 0 0],'LineWidth',2)
plot(main_threshold, volumeEroded,    'Color',[0.7 0.7 0.7],'LineWidth',2)
legend('deposition','net','erosion')

% erosion uncertainty
plot(main_threshold, volumeEroded - uncertaintyVolumeEroded_sum,'--','Color',[0.7 0.7 0.7],'LineWidth',2)
uncertaintyVolumeEroded_sum(volumeEroded + uncertaintyVolumeEroded_sum > 0) = NaN;
plot(main_threshold, volumeEroded + uncertaintyVolumeEroded_sum,'--','Color',[0.7 0.7 0.7],'LineWidth',2)

% deposition uncertainty
plot(main_threshold,volumeDeposited + uncertaintyVolumeDeposited_sum,'--','Color',[0.4 0.4 0.4],'LineWidth',2)
uncertaintyVolumeDeposited_sum(volumeDeposited - uncertaintyVolumeDeposited_sum < 0) = NaN;
plot(main_threshold,volumeDeposited - uncertaintyVolumeDeposited_sum,'--','Color',[0.4 0.4 0.4],'LineWidth',2)

% value at uncertainty level
vmax = max([max(abs(volumeEroded-uncertaintyVolumeEroded_sum)), max(abs(volumeDeposited+uncertaintyVolumeDeposited_sum))]);
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


