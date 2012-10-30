%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function outVal = envmapWorld2Angular(inputEnvMap, dx, dy, dz)
%   Converts an environment map from the world directions to the angular format
%
% Input parameters:
%  - inputEnvMap: light intensities in each input directions
%  - [dx,dy,dz]: world directions
%
% Output parameters:
%  - envMap: environment map in the sky angular format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapWorld2SkyAngular(inputEnvMap, dx, dy, dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Angles

thetaAngular = atan2(dx, dz); % azimuth
phiAngular = atan2(sqrt(dx.^2+dz.^2), dy); % zenith

r = phiAngular./(pi/2);

uAngular = r.*sin(thetaAngular)./2+1/2;
vAngular = 1/2-r.*cos(thetaAngular)./2;

indPos = find(dy>0);

uAngular = uAngular(indPos);
vAngular = vAngular(indPos);

%% Interpolate to get the desired pixel values
envMap = zeros(size(dx,1), size(dx,2), size(inputEnvMap, 3));
for c=1:size(envMap,3)
    o = envMap(:,:,c);
    o(indPos) = interp2(linspace(0,1,size(inputEnvMap,2)), linspace(0,1,size(inputEnvMap,1)), inputEnvMap(:, :, c), uAngular(:), vAngular(:));
    envMap(:, :, c) = o;
end

