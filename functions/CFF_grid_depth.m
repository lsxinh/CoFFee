function [gridEasting,gridNorthing,gridDepth] = CFF_grid_depth(FPBS,field,method,res)
% [gridEasting,gridNorthing,gridDepth] = CFF_grid_depth(FPBS,field,method,res)
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
% - FPBS
% - field: 'DepthZ'
% - method: 'average'
% - res: grid resolution (in m)
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
% 2014-10-13: first version.
%
% EXAMPLE
%
% ...
%
%%%
% Alex Schimel, Deakin University
%%%

% can't do soundings quality method as it uses former datagrams. Might need
% to obtain quality info from new sensors and update this methodology.


% getting data vectors:
Easting = FPBS.Beam.Easting;
Easting = Easting(:);
Northing = FPBS.Beam.Northing;
Northing = Northing(:);
Depth = FPBS.Beam.(field);
Depth = Depth(:);

switch method
    
    case 'average'
        % each sounding has equal weight (method identical to simple
        % average gridding)
        Weight = zeros(size(Depth));
        Weight(~isnan(Depth)) = 1;
        
    case 'angle-weight'
        % beams closest to nadir have higher weight
        Nbeams = max(FPBS.Beam.Rank);
        Weight = ( Nbeams./2 - abs(FPBS.Beam.Rank-Nbeams./2))./Nbeams./2;
        
    case 'soundings-quality'
        % each sounding has a weight associated to the detection method
        % and its quality factors. We define a weighting function that
        % is maximum until the most recurrent value in the data (mode)
        % and then decreases exponentially with soundings quality. We
        % added a -3 factor to decrease the importance of soundings
        % obtained from cross-phase detection in comparison to
        % soundings obtained from amplitude. This methodology and
        % parameters are completely subjective an may be tweaked.
        
        % extract quality data
        Q = File.De_PB_QualityFactor;
        Dec_Pha1_Quality = Q;
        Dec_Pha1_Quality(Q<192) = NaN;
        Dec_Pha1_Quality = Dec_Pha1_Quality-192;
        Dec_Pha2_Quality = Q;
        Dec_Pha2_Quality(Q<128) = NaN;
        Dec_Pha2_Quality(Q>191) = NaN;
        Dec_Pha2_Quality = Dec_Pha2_Quality-128;
        Dec_Amp_NbSamples = Q;
        Dec_Amp_NbSamples(Q>127) = NaN;
        
        % compute weights
        methdiff = 3;
        Dec_Pha1_Weigth = exp(-0.008.*(Dec_Pha1_Quality - mode(Dec_Pha1_Quality(:))+methdiff).^2);
        Dec_Pha1_Weigth(Dec_Pha1_Quality<=mode(Dec_Pha1_Quality(:))-methdiff) = 1;
        Dec_Pha2_Weigth = exp(-0.008.*(Dec_Pha2_Quality - mode(Dec_Pha2_Quality(:))+methdiff).^2);
        Dec_Pha2_Weigth(Dec_Pha2_Quality<=mode(Dec_Pha2_Quality(:))-methdiff) = 1;
        Dec_Amp_Weigth = exp(-0.008.*(Dec_Amp_NbSamples - mode(Dec_Amp_NbSamples(:))).^2);
        Dec_Amp_Weigth(Dec_Amp_NbSamples<=mode(Dec_Amp_NbSamples(:))) = 1;
        
        % combine weights
        Dec_Pha1_Weigth(isnan(Dec_Pha1_Weigth)) = 0;
        Dec_Pha2_Weigth(isnan(Dec_Pha2_Weigth)) = 0;
        Dec_Amp_Weigth(isnan(Dec_Amp_Weigth))   = 0;
        Weight = Dec_Pha1_Weigth + Dec_Pha2_Weigth + Dec_Amp_Weigth;
        
end

% grid
[gridEasting,gridNorthing,gridDepth] = CFF_grid(Easting,Northing,Depth,res,Weight);


