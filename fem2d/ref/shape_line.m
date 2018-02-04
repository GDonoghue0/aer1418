function [shp,shpx] = shape_line(p,x)
% SHAPE_LINE computes the Lagrange polynomials on line segment
% INPUTS:
%   p: polynomial order; must be 1 or 2.
%   x: np by 1 matrix of evaluation points; np is the number of nodes.
% OUTPUTS:
%   shp: np by 2 (or 3) matrix containing Lagrange basis functions
%        evaluated at x.  The number of columns is 2 and 3 for linear and
%        quadratic polynomials, respectively.
%   shpx: np by 2 (or 3) by 2 array containing the derivative of the
%         Lagrange basis functions evaluated at x.  The number of columns 
%         is 3 and 6 for linear and quadratic polynomials, respectively.
%         shp(:,:,1) contains the x-derivative, and shp(:,:,2) constrains
%         the y-derivative.

% Copyright 2018 Masayuki Yano, University of Toronto

xnodes = interp_nodes_line(p);
inv_coeff = monomial_line(p,xnodes);
[psi,psix] = monomial_line(p,x);

shp = psi/inv_coeff;
shpx = psix/inv_coeff;

end


