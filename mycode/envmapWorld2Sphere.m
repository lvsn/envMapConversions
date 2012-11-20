function envMap = envmapWorld2Sphere(inputEnvMap, dx, dy, dz)
% Converts from the world directions to the sphere format.
% 
%	envMap = envmapWorld2Octahedral(inputEnvMap, dx, dy, dz)
%
% Input parameters:
%  - inputEnvMap: light intensities in each input directions
%  - [dx,dy,dz]: world directions
%
% Output parameters:
%  - envMap: environment map in the sky angular format
%
% 
% ----------
% Jean-Francois Lalonde

r = sin(.5.*acos(-dz)) ./ (2.*sqrt(dx.^2+dy.^2));

u = .5+r.*dx;
v = .5-r.*dy;

% Interpolate to get the desired pixel values
dims = size(inputEnvMap, 1);
envMap = zeros(size(dx, 1), size(dx, 2), size(inputEnvMap, 3));
for c=1:size(envMap,3)
    envMap(:, :, c) = reshape(interp2(linspace(0,1,dims), linspace(0,1,dims), inputEnvMap(:, :, c), ...
        u(:), v(:)), size(envMap,1), size(envMap,2));
end
