function [dx,dy,dz,valid] = envmapCube2World(dim)
% Converts from the cube format to the [x,y,z] world directions
%
%   [dx,dy,dz] = envmapCube2World(dim)
%
% Input parameters:
%  - dim: the cube environment dimensions
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
% ----------
% Jean-Francois Lalonde


[u,v] = meshgrid(0:3/(3.*dim-1):3, 0:4/(4.*dim-1):4);
dx = zeros(size(u)); dy = zeros(size(u)); dz = zeros(size(u));

% up
indUp = find(u >= 1 & u < 2 & v < 1);
dx(indUp) = (u(indUp) - 1.5) .* 2;
dy(indUp) = 1;
dz(indUp) = (v(indUp) - 0.5) .* -2;

% left
indLeft = find(u < 1 & v >= 1 & v < 2);
dx(indLeft) = -1;
dy(indLeft) = (v(indLeft) - 1.5) .* -2;
dz(indLeft) = (u(indLeft) - 0.5) .* -2;

% forward
indForward = find(u >= 1 & u < 2 & v >= 1 & v < 2);
dx(indForward) = (u(indForward) - 1.5) .* 2;
dy(indForward) = (v(indForward) - 1.5) .* -2;
dz(indForward) = -1;

% right
indRight = find(u >= 2 & v >= 1 & v < 2);
dx(indRight) = 1;
dy(indRight) = (v(indRight) - 1.5) .* -2;
dz(indRight) = (u(indRight) - 2.5) .* 2;

% down
indDown = find(u >= 1 & u < 2 & v >= 2 & v < 3);
dx(indDown) = (u(indDown) - 1.5) .* 2;
dy(indDown) = -1;
dz(indDown) = (v(indDown) - 2.5) .* 2;

% backward
indBackward = find(u >= 1 & u < 2 & v >= 3);
dx(indBackward) = (u(indBackward) - 1.5) .* 2;
dy(indBackward) = (v(indBackward) - 3.5) .* 2;
dz(indBackward) = 1;

% normalize
norm = sqrt(dx.^2 + dy.^2 + dz.^2);
dx = dx ./ norm; 
dy = dy ./ norm;
dz = dz ./ norm;

if nargout > 3
    valid = ~isnan(dx);
end



