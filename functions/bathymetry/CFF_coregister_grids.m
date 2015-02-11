function [Z1out,Z2out,Xout,Yout] = CFF_coregister_grids(Z1,X1,Y1,Z2,X2,Y2,varargin)
% CFF_new_function(varargin)
%
% DESCRIPTION
%
% This function is to register two grids (X1,Y1,Z1) and (X1,Y1,Z1) on the
% same X,Y by extending the X,Y grids and/or interpolating.
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

% get datasets resolutions:
X1_res = X1(1,2)-X1(1,1);
X2_res = X2(1,2)-X2(1,1);
Y1_res = Y1(1,1)-Y1(2,1);
Y2_res = Y2(1,1)-Y2(2,1);

% output X Y resolution
if nargin==6
    % unspecified, use X1,Y1
    Xout_res = X1_res;
    Yout_res = Y1_res;
elseif nargin==7
    % one value specified, use for both Xout,Yout
    Xout_res = varargin{1};
    Yout_res = varargin{1};
elseif nargin==8
    % two values specified, use first for Xout, second for Yout
    Xout_res = varargin{1};
    Yout_res = varargin{2};
end

% in case the two grids are already co-registered but one file has proper
% easting/northing and not the other.
if all(size(Z1)==size(Z2))
    if X1(1)==1 && Y1(1)==1 && X2(1)~=1 && Y2(1)~=1
        % Z1 has no georeference but is the same size as Z2, assume Z1 has
        % same georeference as Z2.
        X1 = X2;
        Y1 = Y2;
        X1_res = X2_res;
        Y1_res = Y2_res;
    end
    if X2(1)==1 && Y2(1)==1 && X1(1)~=1 && Y1(1)~=1
        % Z2 has no georeference but is the same size as Z1, assume Z2 has
        % same georeference as Z1.
        X2 = X1;
        Y2 = Y1;
        X2_res = X1_res;
        Y2_res = Y1_res;
    end
end

% get coordinates for the max extent
minX = min([X1(1,1);X2(1,1)]);
maxX = max([X1(1,end);X2(1,end)]);
minY = min([Y1(end,1);Y2(end,1)]);
maxY = max([Y1(1,1);Y2(1,1)]);
[Xout,Yout] = meshgrid([minX:Xout_res:maxX],[maxY:-Yout_res:minY]);

% if all resolutions are the same, just fit input grids into output grids
if all(X1_res==[X2_res,Y1_res,Y2_res,Xout_res,Yout_res])
    
    firstcol = find(Xout(1,:)==X1(1));
    firstrow = find(Yout(:,1)==Y1(1));
    Z1out = nan(size(Xout));
    Z1out( firstrow:firstrow+size(Z1,1)-1 , firstcol:firstcol+size(Z1,2)-1 ) = Z1;
    
    firstcol = find(Xout(1,:)==X2(1));
    firstrow = find(Yout(:,1)==Y2(1));
    Z2out = nan(size(Xout));
    Z2out( firstrow:firstrow+size(Z2,1)-1 , firstcol:firstcol+size(Z2,2)-1 ) = Z2;
    
else
    % if not, well we may need some interpolation
    % to code...
end
