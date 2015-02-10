

%% params:
clear all
close all

% files
Z1_file = '.\DATA\WH1_Z_50cm_UTM54S_LAT_p1.tif';
Z2_file = '.\DATA\WH2_Z_50cm_UTM54S_LAT_p1.tif';
Z3_file = '.\DATA\WH3_Z_50cm_UTM54S_LAT_p1.tif';
Z4_file = '.\DATA\WH4_Z_50cm_UTM54S_LAT_p1.tif';

% uncertainty files for LOD
U1_file = '.\DATA\WH1_uncertaintyZ_50cm_UTM54S_p.tif';
U2_file = '.\DATA\WH2_uncertaintyZ_50cm_UTM54S_p.tif';
U3_file = '.\DATA\WH3_uncertaintyZ_50cm_UTM54S_p.tif';
U4_file = '.\DATA\WH4_uncertaintyZ_50cm_UTM54S_p.tif';

% lod method:
lod_method = ['constant',0];    % all data conserved
lod_method = ['constant',0.14]; % constant value
lod_method = ['Ucrit',68];       % t = 1, 68%
lod_method = ['Ucrit',95];       % t = 1.96, 95%

% area of interest
polygon =



%% multi lod analysis
CFF_multi_lod_analysis(Z1_file,Z2_file,polygon,lod_method,U1_file,U2_file);


%% transects:
line = ;

profile = grid_profile({Z1;Z2},line);

