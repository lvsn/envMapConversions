function R = rotationMatrixFromSpinCalcInput(format, input, varargin)
% Helper for converting SpinCalc input to rotation matrices
%
%   e = e.rotate(format, input, <tolerance>, <checks>)
%
%   'format' can be either:
%       - 'DCM'  : direct cosine matrix (input = 3x3 matrix)
%       - 'EA###': euler angles (input = [psi, theta, phi])
%       - 'EV'   : euler rotation vector & angle (input = [vec angle])
%       - 'Q'    : quaternion (input = [q1 q2 q3 q4])
%
% ----------
% Jean-Francois Lalonde

% Get rotation matrix from input
if ~strcmp(format, 'DCM')
    conversion = sprintf('%stoDCM', format);
    R = SpinCalc(conversion, input, varargin{:});
else
    R = input;
end
