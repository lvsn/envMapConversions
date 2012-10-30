%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dx,dy,dz] = envmapSkyAngular2World(dim)
%   Converts an environment map from the sky angular format to the [x,y,z] world directions
%
% Input parameters:
%  - dim: the sky angular environment dimensions
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx,dy,dz] = envmapSkyAngular2World(dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Get the desired world coordinates from the output angular map
[u,v] = meshgrid(0:1/(dim-1):1, 0:1/(dim-1):1);

thetaAngular = atan2(-2.*v+1, 2.*u-1); % azimuth
phiAngular = pi/2.*sqrt((2.*u-1).^2 + (2.*v-1).^2); % zenith

dx = sin(phiAngular).*cos(thetaAngular);
dz = sin(phiAngular).*sin(thetaAngular);
dy = cos(phiAngular);

