function envMap = envmapWorld2SkyAngular(inputEnvMap, dx, dy, dz)
% Converts from the world directions to the "sky angular" format.
% 
%	envMap = envmapWorld2SkyAngular(inputEnvMap, dx, dy, dz)
%
% Input parameters:
%  - inputEnvMap: light intensities in each input directions
%  - [dx,dy,dz]: world directions
%
% Output parameters:
%  - envMap: environment map in the sky angular format
%
% See also:
%   envmapWorld2SkyAngularAngles
% 
% ----------
% Jean-Francois Lalonde

% Get the angles
[uAngular, vAngular, indPos] = envmapWorld2SkyAngularAngles(dx, dy, dz);

% Interpolate to get the desired pixel values
envMap = zeros(size(dx,1), size(dx,2), size(inputEnvMap, 3));
for c=1:size(envMap,3)
    o = envMap(:,:,c);
    o(indPos) = interp2(linspace(0,1,size(inputEnvMap,2)), linspace(0,1,size(inputEnvMap,1)), inputEnvMap(:, :, c), uAngular(:), vAngular(:));
    envMap(:, :, c) = o;
end

