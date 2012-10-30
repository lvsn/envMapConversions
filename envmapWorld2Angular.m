%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function envMap = envmapWorld2Angular(inputEnvMap, dx, dy, dz)
%   Converts an environment map from the world directions to the angular format
%
% Input parameters:
%  - inputEnvMap: light intensities in each input directions
%  - [dx,dy,dz]: world directions
%
% Output parameters:
%  - envMap: environment map in the angular format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapWorld2Angular(inputEnvMap, dx, dy, dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute corresponding image coordinates in the angular map
rAngular = acos(-dz) ./ (2.*pi.*sqrt(dx.^2 + dy.^2));
uAngular = 1/2-rAngular.*dy;
vAngular = 1/2+rAngular.*dx;

%% Interpolate to get the desired pixel values
envMap = zeros(size(dx,1), size(dx,2), size(inputEnvMap, 3));
for c=1:size(envMap,3)
    envMap(:, :, c) = reshape(interp2(0:1/(size(inputEnvMap,2)-1):1, 0:1/(size(inputEnvMap,1)-1):1, inputEnvMap(:, :, c)', uAngular(:), vAngular(:)), size(envMap,1), size(envMap,2));
end

