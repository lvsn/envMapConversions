%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function envMap = envmapWorld2Cube(inputEnvMap, dx, dy, dz)
%   Converts an environment map from the world directions to the cube map format
%
% Input parameters:
%  - inputEnvMap: light intensities in each input directions
%  - [dx,dy,dz]: world directions
%
% Output parameters:
%  - envMap: environment map in the cube format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function envMap = envmapWorld2Cube(inputEnvMap, dx, dy, dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2006-2009 Jean-Francois Lalonde
% Carnegie Mellon University
% Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get u-v coordinates
uCube = zeros(size(dx));
vCube = zeros(size(dx));

% forward
indForward = find(dz <= 0 & dz <= -abs(dx) & dz <= -abs(dy));
uCube(indForward) = 1.5 - 0.5 .* dx(indForward) ./ dz(indForward);
vCube(indForward) = 1.5 + 0.5 .* dy(indForward) ./ dz(indForward);

% backward
indBackward = find(dz >= 0 & dz >= abs(dx) & dz >= abs(dy));
uCube(indBackward) = 1.5 + 0.5 .* dx(indBackward) ./ dz(indBackward);
vCube(indBackward) = 3.5 + 0.5 .* dy(indBackward) ./ dz(indBackward);

% down
indDown = find(dy <= 0 & dy <= -abs(dx) & dy <= -abs(dz));
uCube(indDown) = 1.5 - 0.5 .* dx(indDown) ./ dy(indDown);
vCube(indDown) = 2.5 - 0.5 .* dz(indDown) ./ dy(indDown);

% up
indUp = find(dy >= 0 & dy >= abs(dx) & dy >= abs(dz));
uCube(indUp) = 1.5 + 0.5 .* dx(indUp) ./ dy(indUp);
vCube(indUp) = 0.5 - 0.5 .* dz(indUp) ./ dy(indUp);

% left
indLeft = find(dx <= 0 & dx <= -abs(dy) & dx <= -abs(dz));
uCube(indLeft) = 0.5 + 0.5 .* dz(indLeft) ./ dx(indLeft);
vCube(indLeft) = 1.5 + 0.5 .* dy(indLeft) ./ dx(indLeft);

% right
indRight = find(dx >= 0 & dx >= abs(dy) & dx >= abs(dz));
uCube(indRight) = 2.5 + 0.5 .* dz(indRight) ./ dx(indRight);
vCube(indRight) = 1.5 - 0.5 .* dy(indRight) ./ dx(indRight);

%% Interpolate to get the desired pixel values
envMap = zeros(size(dx,1), size(dx,2), size(inputEnvMap, 3));
for c=1:3
    envMap(:, :, c) = reshape(interp2(0:3/(size(inputEnvMap,2)-1):3, 0:4/(size(inputEnvMap,1)-1):4, inputEnvMap(:, :, c), uCube(:), vCube(:)), size(envMap,1), size(envMap,2));
end


