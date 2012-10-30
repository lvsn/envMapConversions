%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dx,dy,dz] = envmapAngular2World(dim)
%   Converts an environment map from the angular format to the [x,y,z] world directions
%
% Input parameters:
%  - dim: the angular environment dimensions
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx,dy,dz] = envmapAngular2World(dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the desired world coordinates from the output angular map
[u,v] = meshgrid(0:1/(dim-1):1, 0:1/(dim-1):1);

thetaAngular = atan2(-2.*v+1, 2.*u-1);
phiAngular = pi.*sqrt((2.*u-1).^2 + (2.*v-1).^2);

dx = sin(phiAngular).*cos(thetaAngular);
dy = sin(phiAngular).*sin(thetaAngular);
dz = -cos(phiAngular);