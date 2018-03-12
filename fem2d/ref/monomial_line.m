function [shp,shpx,shpxx] = monomial_line(p,x)
% MONOMIAL_LINE computes monomials in 1d
% INPUTS:
%   p: polynomial order; must be 1 or 2.
%   x: np by 1 matrix of evaluation points; np is the number of nodes.
% OUTPUTS:
%   shp: np by 2 (or 3) matrix containing monomial basis functions
%        evaluated at x.  The number of columns is 2 and 3 for linear and
%        quadratic polynomials, respectively.
%   shpx: np by 2 (or 3) matrix containing the derivative of the
%         monomial basis functions evaluated at x.  The number of columns 
%         is 2 and 3 for linear and quadratic polynomials, respectively.
%   shpxx: np by 2 (or 3) matrix containing the second derivative of the
%          monomial basis functions evlauated at x.

% Copyright 2018 Masayuki Yano, University of Toronto

np = size(x,1);
oo = ones(np,1);
zz = zeros(np,1);
switch p
    case 1
        shp = [oo, x];
        shpx = [zz, oo];
        shpxx = [zz, zz];
    case 2
        shp = [oo, x, x.^2];
        shpx = [zz, oo, 2.*x];
        shpxx = [zz, zz, 2.*oo];
    otherwise
        error('unsupported polynomial order');
end
end