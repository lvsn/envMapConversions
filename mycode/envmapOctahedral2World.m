function [dx,dy,dz,valid] = envmapOctahedral2World(dim)
% Converts from the octahedral format to the [x,y,z] world directions
% 
%   [dx,dy,dz,valid] = envmapOctahedral2World(dim)
%   
%
% Input parameters:
%  - dim: the sky angular environment dimensions
%
% Output parameters:
%  - [dx,dy,dz]: world directions
%
% ----------
% Jean-Francois Lalonde

% Get the desired world coordinates from the output angular map
[u,v] = meshgrid(linspace(-1,1,dim), linspace(-1,1,dim));

thetaAngular = zeros(size(u));
phiAngular = zeros(size(u));

% split into triangles
% theta = zenith
% phi = azimuth

% center part (upper hemisphere)
t1 = u>=0 & v>=0 & v<1-u;
t2 = u<0 & v>=0 & v<u+1;
t3 = u<0 & v<0 & v>=-u-1;
t4 = u>=0 & v<0 & v>=u-1;

% outer part (lower hemisphere)
t5 = u>=0 & v>=0 & v>=1-u;
t6 = u<0 & v>= 0 & v>=u+1;
t7 = u<0 & v<0 & v<-u-1;
t8 = u>=0 & v<0 & v<u-1;

thetaAngular(t1|t5) = (v(t1|t5)+u(t1|t5))*pi/2;
thetaAngular(t2|t6) = (v(t2|t6)-u(t2|t6))*pi/2;
thetaAngular(t3|t7) = (-u(t3|t7)-v(t3|t7))*pi/2;
thetaAngular(t4|t8) = (u(t4|t8)-v(t4|t8))*pi/2;

phiAngular(t1|t2|t3|t4) = atan2(v(t1|t2|t3|t4),u(t1|t2|t3|t4));
phiAngular(t5) = atan2(1-u(t5), 1-v(t5));
phiAngular(t6) = atan2(u(t6)+1, v(t6)-1);
phiAngular(t7) = -atan2(u(t7)+1, -1-v(t7));
phiAngular(t8) = -atan2(1-u(t8), 1+v(t8));

dx = cos(pi/2-thetaAngular).*cos(phiAngular);
dz = cos(pi/2-thetaAngular).*sin(phiAngular);
dy = sin(pi/2-thetaAngular);


% dx = sin(thetaAngular).*cos(phiAngular);
% dz = sin(thetaAngular).*sin(phiAngular);
% dy = cos(thetaAngular);

if nargout > 3
    % everything is valid
    valid = true(size(u));
end

