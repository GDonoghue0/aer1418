function [shp,shpx] = monomial(p,x)
% MONOMIAL computes the (complete) monomials in 2d
% INPUTS:
%   p: polynomial order; must be 1 or 2.
%   x: np by 2 matrix of evaluation points; np is the number of nodes.
% OUTPUTS:
%   shp: np by 3 (or 6) matrix containing monomial basis functions
%        evaluated at x.  The number of columns is 3 and 6 for linear and
%        quadratic polynomials, respectively.
%   shpx: np by 3 (or 6) by 2 array containing the derivative of the
%         monomial basis functions evaluated at x.  The number of columns 
%         is 3 and 6 for linear and quadratic polynomials, respectively.
%         shp(:,:,1) contains the x-derivative, and shp(:,:,2) constrains
%         the y-derivative.
np = size(x,1);
oo = ones(np,1);
zz = zeros(np,1);
switch p
    case 1
        shp = [oo, x(:,1), x(:,2)];
        shpx = zeros(np,3,2);
        shpx(:,:,1) = [zz, oo, zz];
        shpx(:,:,2) = [zz, zz, oo];
    case 2
        shp = [oo, x(:,1), x(:,2), x(:,1).^2, x(:,1).*x(:,2), x(:,2).^2];
        shpx = zeros(np,6,2);
        shpx(:,:,1) = [zz, oo, zz, 2.*x(:,1), x(:,2), zz];
        shpx(:,:,2) = [zz, zz, oo, zz, x(:,1), 2.*x(:,2)];
    otherwise
        error('unsupported polynomial order');
end
end