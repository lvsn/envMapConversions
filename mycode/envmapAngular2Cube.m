%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function envMap = envmapLatLong2ACube(angularEnvMap, dim)
%   Converts an environment map from the angular format to the "cube" format.
%
% Input parameters:
%  - angularEnvMap: environment map in angular format
%  - dim: dimensions of the output environment map (dim x dim)
%
% Output parameters:
%  - envMap: environment map in cube format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapAngular2Cube(angularEnvMap, dim)
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
envMap(indInside,:) = envmapWorld2Angular(angularEnvMap, dxCube(indInside), dyCube(indInside), dzCube(indInside));
envMap = reshape(envMap, size(dxCube,1), size(dxCube,2), 3);


