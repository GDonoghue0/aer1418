function A = circle_area(p,h)
% CIRCLE_AREA is a driver file to test area integration
%   The driver file demonstrates the following idea: (i) isoparametric
%   mapping; (ii) convergence of geometry approximation.
% INPUTS:
%   p: polynomial degree (must be 1 or 2)
%   h: characteristic size of elements
% OUTPUT:
%   A: area of the domain

% Copyright 2018 Masayuki Yano, University of Toronto

% discretization parameters
%pquad = ?;

% make reference
%ref = make_ref_tri(p,pquad);

% generate mesh
%mesh = make_circle_mesh(h, p==2); 

% useful parameters
ntri = size(mesh.tri,1);

% sum area
A = 0.0;
for elem = 1:ntri
    % TODO: compute the determinant of the Jacobian 
    
    % TODO: integrate det(J) and add to A    
end

end