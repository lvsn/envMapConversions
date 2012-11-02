function [dx,dy,dz,valid] = envmapLatLong2World(height)
% Converts from the angular format to the [x,y,z] world directions
%
%   [dx,dy,dz,valid] = envmapLatLong2World(height)
%
% Input parameters:
%  - height: height of the latitude-longitude output (width = 2*height)
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
% ----------
% Jean-Francois Lalonde

%% Get the desired world coordinates from the output lat-long map
[u,v] = meshgrid(0:2/(2.*height-1):2, 0:1/(height-1):1);

thetaLatLong = pi.*(u-1);
phiLatLong = pi.*v;

dx = sin(phiLatLong).*sin(thetaLatLong);
dy = cos(phiLatLong);
dz = -sin(phiLatLong).*cos(thetaLatLong);

% all the coordinates are valid
if nargout > 3
    valid = true(size(dx));
end
