%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function envMap = envmapSkyAngular2LatLong(angularEnvMap, dim)
%   Converts an environment map from the "sky angular" format to the lat-long format
%
% Input parameters:
%  - angularEnvMap: environment map in sky angular format
%  - dim: dimensions of the output environment map (dim x 2*dim)
%
% Output parameters:
%  - envMap: environment map in the lat-long format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapSkyAngular2LatLong(angularEnvMap, dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the 3-D directions
[dxLatLong, dyLatLong, dzLatLong] = envmapLatLong2World(dim);

%% Get the angular environment map representation
envMap = envmapWorld2SkyAngular(angularEnvMap, dxLatLong, dyLatLong, dzLatLong);
