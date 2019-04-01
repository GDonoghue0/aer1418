function [shp,shpx] = monomial_tri(p,x)
% MONOMIAL_TRI computes the (complete) monomials in 2d
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

% Copyright 2018 Masayuki Yano, University of Toronto

end
