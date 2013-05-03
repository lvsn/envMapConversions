function omega = tetrahedronSolidAngle(a, b, c)
% Computes the solid angle subtended by a tetrahedron.
%
%   omega = tetrahedronSolidAngle(a, b, c)
%
% The tetrahedron is defined by three vectors (a, b, c) which define the
% vertices of the triangle with respect to an origin.
%
% For more details, see:
%   http://en.wikipedia.org/wiki/Solid_angle#Tetrahedron  
%
% We use L'Huilier's theorem because it's simpler to parallelize
% 
% ----------
% Jean-Francois Lalonde

assert(size(a, 1) == 3, 'a must be a 3xN matrix');
assert(size(b, 1) == 3, 'b must be a 3xN matrix');
assert(size(c, 1) == 3, 'c must be a 3xN matrix');

theta_a = acos(sum(b.*c, 1));
theta_b = acos(sum(a.*c, 1));
theta_c = acos(sum(a.*b, 1));

theta_s = (theta_a + theta_b + theta_c)./2;

product = tan(theta_s/2).*tan((theta_s-theta_a)./2).*...
    tan((theta_s-theta_b)./2).*tan((theta_s-theta_c)./2);

% account for negatives (not sure what to do here?)
product(product<0) = 0;
tmp = sqrt(product);

omega = 4.*atan(tmp);

return;


a = column(a(:,51));
b = column(b(:,51));
c = column(c(:,51));

al = norm(a); a = a./al;
bl = norm(b); b = b./bl;
cl = norm(c); c = c./cl;

determ = det([a b c]);

div = al*bl*cl + (a'*b)*cl + (a'*c)*bl + (b'*c)*al;
at = atan2(determ, div);

at(at<0) = at(at<0) + pi; % If det > 0 and div < 0 arctan2 returns < 0, so add pi.

omega = 2*at;
    
 
