%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function envMap = envmapCube2Angular(cubeEnvMap, dim)
%   Converts an environment map from the cube format to the angular format.
%
% Input parameters:
%  - latLongEnvMap: environment map in latitude-longitude format
%  - dim: dimensions of the output environment map (dim x dim)
%
% Output parameters:
%  - envMap: environment map in angular format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapCube2Angular(cubeEnvMap, dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the 3-D directions
[dxAngular, dyAngular, dzAngular] = envmapAngular2World(dim);

%% Get the values from the lat-long environment map representation
envMap = envmapWorld2Cube(cubeEnvMap, dxAngular, dyAngular, dzAngular);

%% Mask out the outside of the circle
circleMask = zeros(dim);
[c,r] = meshgrid(1:dim,1:dim);

indCircle = (c-dim/2).^2+(r-dim/2).^2 <= (dim/2).^2;
circleMask(indCircle) = 1;

envMap = envMap .* repmat(circleMask, [1 1 3]);