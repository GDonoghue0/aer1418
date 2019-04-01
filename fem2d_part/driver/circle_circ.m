function C = circle_circ(p,h)
% CIRCLE_CIRC is a driver file to test boundary integration
%   The driver file demonstrates the following idea: (i) isoparametric
%   mapping; (ii) convergence of geometry approximation.

% Copyright 2018 Masayuki Yano, University of Toronto

% discretization parameters
% pquad = ?;

% make reference
% ref = make_ref_tri(p,pquad);

% generate mesh and create boundary group
% mesh = make_circle_mesh(h, p==2);
% mesh = make_bgrp(mesh);

% sum circumference
C = 0.0;
for bgrp = 1:length(mesh.bgrp)
    for edge = 1:size(mesh.bgrp{bgrp},1)
        % get element, local edge, local nodes, and global nodes
        elem = mesh.bgrp{bgrp}(edge,3);
        ledge = mesh.bgrp{bgrp}(edge,4);
        lnode = ref.f2n(:,ledge);
        tril = mesh.tri(elem,lnode).';
        
        % TODO: compute the determinant of the facet Jacobian
        
        % TODO: integrate det(J) and add to C
    end
end

end
