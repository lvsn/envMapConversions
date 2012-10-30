%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function outVal = envmapWorld2LatLong(inputEnvMap, dx, dy, dz)
%   Converts an environment map from the world directions to the lat-long directions
%
% Input parameters:
%  - inputEnvMap: light intensities in each input directions
%  - [dx,dy,dz]: world directions
%
% Output parameters:
%  - envMap: environment map value at the [x,y,z] directions specified
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapWorld2LatLong(inputEnvMap, dx, dy, dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute corresponding image coordinates in the lat-long map
uLatLong = 1 + (1/pi) .* atan2(dx, -dz);
vLatLong = (1/pi) .* acos(dy);

%% Interpolate to get the desired pixel values
envMap = zeros(size(dx, 1), size(dx, 2), size(inputEnvMap, 3));
for c=1:size(envMap,3)
    envMap(:, :, c) = reshape(interp2(0:2/(size(inputEnvMap,2)-1):2, 0:1/(size(inputEnvMap,1)-1):1, inputEnvMap(:, :, c), uLatLong(:), vLatLong(:)), size(envMap,1), size(envMap,2));
end

