%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function envMap = envmapLatLong2Cube(latLongEnvMap, dim)
%   Converts an environment map from the lat-long format to the cube format.
%
% Input parameters:
%  - latLongEnvMap: environment map in latitude-longitude format
%  - dim: dimensions of the output environment map (dim x dim)
%
% Output parameters:
%  - envMap: environment map in cube format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapLatLong2Cube(latLongEnvMap, dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the 3-D directions
[dxCube, dyCube, dzCube] = envmapCube2World(dim);
indInside = (~isnan(dxCube) & ~isnan(dyCube) & ~isnan(dzCube));

%% Get the values from the lat-long environment map representation
envMap = zeros(size(dxCube, 1)*size(dxCube, 2), 3);
envMap(indInside,:) = envmapWorld2LatLong(latLongEnvMap, dxCube(indInside), dyCube(indInside), dzCube(indInside));
envMap = reshape(envMap, size(dxCube,1), size(dxCube,2), 3);


