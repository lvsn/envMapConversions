%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dx,dy,dz] = envmapLatLong2World(height)
%   Converts an environment map from the angular format to the [x,y,z] world directions
%
% Input parameters:
%  - height: height of the latitude-longitude output (width = 2*height)
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx,dy,dz] = envmapLatLong2World(height)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the desired world coordinates from the output lat-long map
[u,v] = meshgrid(0:2/(2.*height-1):2, 0:1/(height-1):1);

thetaLatLong = pi.*(u-1);
phiLatLong = pi.*v;

dx = sin(phiLatLong).*sin(thetaLatLong);
dy = cos(phiLatLong);
dz = -sin(phiLatLong).*cos(thetaLatLong);