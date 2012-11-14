function envMap = envmapWorld2Octahedral(inputEnvMap, dx, dy, dz)
% Converts from the world directions to the octahedral format.
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
% See also:
%   envmapWorld2SkyAngularAngles
% 
% ----------
% Jean-Francois Lalonde

% Map (x,y,z) to (u,v) in the octahedral map
phi = atan2(dz, dx); % azimuth
theta = atan2(sqrt(dx.^2+dz.^2), dy); % zenith

uOctahedral = zeros(size(phi));
vOctahedral = zeros(size(phi));

% center part (upper hemisphere)
t1 = theta<pi/2 & phi>=0 & phi<pi/2;
t2 = theta<pi/2 & phi>=pi/2;
t3 = theta<pi/2 & phi<-pi/2;
t4 = theta<pi/2 & phi<0 & phi>=-pi/2;

% outer part (lower hemisphere)
t5 = theta>=pi/2 & phi>=0 & phi<pi/2;
t6 = theta>=pi/2 & phi>=pi/2;
t7 = theta>=pi/2 & phi<-pi/2;
t8 = theta>=pi/2 & phi<0 & phi>=-pi/2;

thetap = 2*theta/pi;

uOctahedral(t1) = thetap(t1)./(1+tan(phi(t1)));
uOctahedral(t2) = thetap(t2)./(tan(phi(t2))-1);
uOctahedral(t3) = thetap(t3)./(-tan(phi(t3))-1);
uOctahedral(t4) = thetap(t4)./(1-tan(phi(t4)));

uOctahedral(t5) = ((thetap(t5)-1).*tan(phi(t5))+1)./(tan(phi(t5))+1);
uOctahedral(t6) = ((thetap(t6)-1).*tan(phi(t6))-1)./(1-tan(phi(t6)));
uOctahedral(t7) = -((thetap(t7)-1).*tan(phi(t7))+1)./(tan(phi(t7))+1);
uOctahedral(t8) = ((thetap(t8)-1).*tan(phi(t8))-1)./(tan(phi(t8))-1);

vOctahedral(t1|t5) = thetap(t1|t5)-uOctahedral(t1|t5);
vOctahedral(t2|t6) = thetap(t2|t6)+uOctahedral(t2|t6);
vOctahedral(t3|t7) = -thetap(t3|t7)-uOctahedral(t3|t7);
vOctahedral(t4|t8) = -thetap(t4|t8)+uOctahedral(t4|t8);


% phiAngular(t5) = atan2(1-u(t5), 1-v(t5));
% phiAngular(t6) = atan2(u(t6)+1, v(t6)-1);
% phiAngular(t7) = -atan2(u(t7)+1, -1-v(t7));
% phiAngular(t8) = -atan2(1-u(t8), 1+v(t8));

% thetaAngular(t1|t5) = (v(t1|t5)+u(t1|t5))*pi/2;
% thetaAngular(t2|t6) = (v(t2|t6)-u(t2|t6))*pi/2;
% thetaAngular(t3|t7) = (-u(t3|t7)-v(t3|t7))*pi/2;
% thetaAngular(t4|t8) = (u(t4|t8)-v(t4|t8))*pi/2;

% Interpolate to get the desired pixel values
dims = size(inputEnvMap, 1);
envMap = zeros(size(dx, 1), size(dx, 2), size(inputEnvMap, 3));
for c=1:size(envMap,3)
    envMap(:, :, c) = reshape(interp2(linspace(-1,1,dims), linspace(-1,1,dims), inputEnvMap(:, :, c), ...
        uOctahedral(:), vOctahedral(:)), size(envMap,1), size(envMap,2));
end

