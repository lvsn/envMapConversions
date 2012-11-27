function [dx,dy,dz,valid] = envmapSphere2World(dim)
% Converts from the sphere format to the [x,y,z] world directions
% 
%   [dx,dy,dz,valid] = envmapSphere2World(dim)
%   
%
% Input parameters:
%  - dim: the sphere environment dimensions
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
% Note: see p. 402 from HDRI book. 
% There's a typo in the book: dx is flipped with dy.
%
% ----------
% Jean-Francois Lalonde

% Get the desired world coordinates from the output angular map
[u,v] = meshgrid(linspace(-1,1,dim), linspace(-1,1,dim));

r = sqrt(u.^2 + v.^2);
theta = atan2(u, -v);

phi = zeros(size(theta));
valid = r<=1;
phi(valid) = 2.*asin(r(valid));

dx = sin(phi).*sin(theta);
dy = sin(phi).*cos(theta);
dz = -cos(phi);

