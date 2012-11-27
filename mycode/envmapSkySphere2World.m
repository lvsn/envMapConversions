function [dx,dy,dz,valid] = envmapSkySphere2World(dim)
% Converts from the "sky-only" sphere format to the [x,y,z] world directions
% 
%   [dx,dy,dz,valid] = envmapSkySphere2World(dim)
%   
%
% Input parameters:
%  - dim: the sphere environment dimensions
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
% ----------
% Jean-Francois Lalonde

% Get the desired world coordinates from the output angular map
[u,v] = meshgrid(linspace(-1,1,dim), linspace(-1,1,dim));
valid = u.^2 + v.^2 <= 1;

uM = u*sqrt(2)/2;
vM = v*sqrt(2)/2;

r = sqrt(uM.^2 + vM.^2);

phi = zeros(size(valid));
theta = zeros(size(valid));

theta(valid) = atan2(uM(valid), -vM(valid));
phi(valid) = 2.*asin(r(valid));

% dx = sin(phi).*sin(theta);
% dy = sin(phi).*cos(theta);
% dz = -cos(phi);

dx = sin(phi).*sin(theta);
dy = cos(phi);
dz = sin(phi).*cos(theta);


