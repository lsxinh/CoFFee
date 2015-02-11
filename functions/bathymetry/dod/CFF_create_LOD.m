function [LOD,X,Y] = CFF_create_LOD(varargin)
% [LOD,LOD_easting,LOD_northing] = CFF_create_LOD(varargin)
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

% methods:

% lod_method = ['constant',0   ,X,Y]; % all data conserved
% lod_method = ['constant',0.14,X,Y]; % constant value
% lod_method = ['Ucrit'   ,68  ,U1,X1,Y1,U2,X2,Y2];       % t = 1, 68%
% lod_method = ['Ucrit'   ,95  ,U1,X1,Y1,U2,X2,Y2];       % t = 1.96, 95%

method = varargin{1};

switch method
    case 'constant'
        
        threshold = varargin{2};
        
        
            
    case 'Ucrit'
        
        % read rest of parameters
        
        CONF = varargin{2};
        
        U1 = varargin{3};
        X1 = varargin{4};
        Y1 = varargin{5};
        
        U2 = varargin{6};
        X2 = varargin{7};
        Y2 = varargin{8};
        
        % coregister grids
        [U1,U2,X,Y] = CFF_coregister_grids(U1,X1,Y1,U2,X2,Y2);
        
        % t value corresponding to an confidence limit in percentage:
        t = CFF_critical_z_value(CONF);
        
        % Uncertainty propagated in quadrature:
        UD = sqrt( U1.^2 + U2.^2 );
        
        % Ucrit (as per Brasington et al., 2003; Lane et al., 2003; Wheaton et al., 2010; Milan et al., 2011; Lallias-Tacon et al. 2013; etc.
        % ="minimum level of detection threshold (LOD)"
        % ="LoD grid" (Carley et al., 2012)
        LOD = t.*UD;
        
end




