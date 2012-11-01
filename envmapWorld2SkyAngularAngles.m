function [uSkyAngular, vSkyAngular, indPos] = envmapWorld2SkyAngularAngles(dx, dy, dz)
% Converts (x,y,z) coordinates to "sky angular" format
%
%   [uSkyAngular, vSkyAngular] = envmapWorld2SkyAngularAngles(dx, dy, dz)
%
%   The inputs 'dx', 'dy', 'dz' are (x,y,z) coordinates in the world 
%   reference frame.
%
%   The outputs 'uSkyAngular' and 'vSkyAngular' are the (u,v) image 
%   coordinates, in the [0,1] interval.
% 
% See also:
%   envmapWorld2SkyAngular
% 
% ----------
% Jean-Francois Lalonde

thetaAngular = atan2(dx, dz); % azimuth
phiAngular = atan2(sqrt(dx.^2+dz.^2), dy); % zenith

r = phiAngular./(pi/2);

uSkyAngular = r.*sin(thetaAngular)./2+1/2;
vSkyAngular = 1/2-r.*cos(thetaAngular)./2;

indPos = find(dy>0);

uSkyAngular = uSkyAngular(indPos);
vSkyAngular = vSkyAngular(indPos);

